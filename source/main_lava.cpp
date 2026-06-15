#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>

#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>
#include <fstream>

#include "json.hpp"
#include "Taylor_Galerkin_lava.h"

using json = nlohmann::json;
using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = distributed_vector; // Typedef for local q_0 vector

static constexpr char VARNAME_1[255] = "dem"; 

static std::vector<double> dem;

double h_min, L, H, res;
int NUM_REFINEMENTS, Nx, Ny;


// Refinement rule
static int uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int hanging_refinement (tmesh::quadrant_iterator quadrant)
{
    double x_minus, x_plus, y_minus, y_plus;
    x_minus = quadrant->p (0, 0);
    x_plus  = quadrant->p (0, 1);
    y_minus = quadrant->p (1, 0);
    y_plus  = quadrant->p (1, 2);

    double x_center, y_center;
    x_center = quadrant->centroid (0);
    y_center = quadrant->centroid (1);

    const auto marker = x_center<3./4.*L && x_center>L/4. && y_center<3./4.*H && y_center>H/4.;
    return marker; 
}


inline static int raster_2_vector(const int i_x, const int i_y)
{ return(i_y + Ny*i_x); }

static std::array<int,3> global_coord_2_raster(const double x,
        const double y)
{
    // nearest neighbor
    const int i_x = static_cast<int>(std::round(x / res));
    const int i_y = static_cast<int>(-std::round(y / res)) + (Ny-1);
    const int ii  = raster_2_vector(i_x, i_y);

    return(std::array<int,3>{{ ii, i_x, i_y }});
}

inline static double raster_value(const double x, const double y)
{
    // bilinear interp.
    const double ix = x/res;
    const double iy = -y/res + (Ny-1);

    const std::array<int,2> ix_q = {static_cast<int>(std::floor(ix)), 
        static_cast<int>(std::ceil (ix))}; 
    const std::array<int,2> iy_q = {static_cast<int>(std::floor(iy)), 
        static_cast<int>(std::ceil (iy))};

    const auto Dx_adi = ix_q[1]-ix_q[0];
    const auto Dy_adi = iy_q[1]-iy_q[0];

    std::array<double,4> gamma = {0,0,0,0};
    for (int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++)
        {
            gamma[i+j*2] = dem[raster_2_vector(ix_q[(i+1)%2],iy_q[(j+1)%2])];

            gamma[i+j*2] *= Dx_adi!= 0 ? std::abs(ix_q[i]-ix)/Dx_adi : .5;
            gamma[i+j*2] *= Dy_adi!= 0 ? std::abs(iy_q[j]-iy)/Dy_adi : .5;
        }
    }

    return(gamma[0]+gamma[1]+gamma[2]+gamma[3]);
}

inline double dem_fun (const double xx, const double yy)
{
#if SET_TEST == 11 // WB test for the set of nonlinear equations
    return (xx>3 && xx<7 && yy>3 && yy<7) ? 5.0 * std::exp(-2.0/5.0 * (std::pow(xx-L/2.0,2) + std::pow(yy-H/2.0,2))) : 0;
#endif
    return 0;
    //return(raster_value(xx,yy));
}

#if SET_TEST == 12 
// Parameters
static constexpr double h0     = 1.0;
static constexpr double hmin   = 0.9;
static constexpr double x_0    = 0.5;
static constexpr double y_0    = 0.5;
static constexpr double r0     = 0.25;
static constexpr double u_inf  = 6.0;

auto Gamma = [](double g) {
    return M_PI/2./r0*std::sqrt((g*(h0-hmin))/(3.*M_PI*M_PI/64.-1./4.));
};
 
auto omega = [](double r, double g, double rho) -> double {
    return 2*Gamma(g)*std::cos(rho/2.)*std::cos(rho/2.)*(r<=r0);
};
#endif


inline double h0_fun (const double xx, const double yy, const double g) 
{
    double h_ini = 0;
    
#if SET_TEST == 11
    h_ini = 10. - dem_fun(xx,yy);
#elif SET_TEST == 12 
    double r     = std::sqrt((xx-x_0)*(xx-x_0) + (yy-y_0)*(yy-y_0));
    double rho   = M_PI*r/r0; 
    double corr = 4.;
    auto H_vortex = [](double s) { 
        return std::cos(2.*s)/8. + s*std::sin(2.*s)/4. + std::cos(2.*s)*std::cos(2.*s)/64. + 3.*s*s/16. + s*std::cos(2.*s)*std::sin(2.*s)/16.;
    };
    h_ini = h0 - (r<=r0)*1./g*(2.*Gamma(g)*r0/M_PI)*(2.*Gamma(g)*r0/M_PI)*(H_vortex(M_PI/2.) - H_vortex(rho/2.))*corr;
#elif SET_TEST == 21 
    h_ini = 1. + 3.*std::exp(-5.*( (xx-L/4.)*(xx-L/4.) )/(.1*L)/(.1*L));
#endif

    return (h_ini);
} 

inline double Ux0_fun (const double xx, const double yy, const double g)  
{ 
    double U_ini = 0.;

#if SET_TEST == 12
    double r       = std::sqrt((xx-x_0)*(xx-x_0) + (yy-y_0)*(yy-y_0));
    double rho     = M_PI*r/r0;
    U_ini   = h0_fun(xx,yy,g)*(u_inf - (yy - y_0)*omega(r,g,rho));
#endif

    return (U_ini);
}

inline double Uy0_fun (const double xx, const double yy, const double g) 
{ 
    double V_ini = 0.;

#if SET_TEST == 12
    double r       = std::sqrt((xx-x_0)*(xx-x_0) + (yy-y_0)*(yy-y_0));
    double rho     = M_PI*r/r0;
    V_ini   = h0_fun(xx,yy,g)*(xx - x_0)*omega(r,g,rho);
#endif

    return (V_ini);
}

inline double Th0_fun (double xx, double yy)
{
#if SET_TEST == 11
    return 1; // in km \cdot K
#endif
    return 0.;
}


template <class T> 
void quadrant_marker_list (tmesh::quadrant_iterator& q, 
        const T& only_h, const T& only_Ux, const T& only_Uy, const double& dt, std::set<int>& output)
{
    std::array<double,4> h_current = {0,0,0,0};
    std::array<double,2> vel       = {0,0};
    for (int ii = 0; ii < 4; ++ii)
    {
        if (! q->is_hanging (ii)){
            const auto & h_candidate = only_h[q->gt (ii)];
            h_current[ii] = h_candidate;

            vel[0] += h_candidate>h_min ? only_Ux[q->gt (ii)]/h_candidate : 0.;
            vel[1] += h_candidate>h_min ? only_Uy[q->gt (ii)]/h_candidate : 0.;
        }
        else
        {
            const auto & h_candidate = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
            h_current[ii] = h_candidate;

            vel[0] += h_candidate>h_min ? .5*(only_Ux[q->gparent (0,ii)] + only_Ux[q->gparent (1,ii)])/h_candidate : 0.;

            vel[1] += h_candidate>h_min ? .5*(only_Uy[q->gparent (0,ii)] + only_Uy[q->gparent (1,ii)])/h_candidate : 0.;
        }
    }
    vel[0] /= 4.;
    vel[1] /= 4.;

    const bool basin_check = ((h_current[0]+h_current[1]+h_current[2]+h_current[3])> h_min && 
            (h_current[0]*h_current[1]*h_current[2]*h_current[3])<=h_min) ? true : false;

    if (basin_check)
    {
        std::vector<std::tuple<tmesh::quadrant_iterator, int> > quadrant_list;
        output.insert(q->get_global_quad_idx ());
        for (auto quadrant_nei  = q->begin_neighbor_sweep();
                quadrant_nei != q->end_neighbor_sweep (); ++quadrant_nei)
        {
            quadrant_list.push_back(std::tuple<tmesh::quadrant_iterator, int>{quadrant_nei, quadrant_nei->get_global_quad_idx ()});
        }

        const double xx_ini = q->centroid (0);
        const double yy_ini = q->centroid (1);

        const double xx_fin = xx_ini + vel[0]*dt;
        const double yy_fin = yy_ini + vel[1]*dt;

        for (auto & qq : quadrant_list)
        {
            auto & quadrant = std::get<0>(qq);

            const double x1 = quadrant->p(0,0);
            const double x2 = quadrant->p(0,1);
            const double y1 = quadrant->p(1,0);
            const double y2 = quadrant->p(1,2); 

            const bool cond1 = ( yy_fin - y1)>=0;
            const bool cond2 = (-yy_fin + y2)>=0;
            const bool cond3 = ( xx_fin - x1)>=0;
            const bool cond4 = (-xx_fin + x2)>=0;

            if (cond1 && cond2 && cond3 && cond4) // the point is internal to the current quadrant
            {
                output.insert(std::get<1>(qq));
            }
        }
    }
}


// Re-Define tic and toc to add an MPI_Barrier
#define TIC()  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int main (int argc, char **argv)
{
    // parse input file,
    MPI_Init (&argc, &argv);
    std::ifstream input_file(argv[1]);
    json input_data = json::parse(input_file);

    res                                                     = input_data["raster resolution"];
    Nx                                                      = input_data["number raster columns"];
    Ny                                                      = input_data["number raster rows"];
    NUM_REFINEMENTS                                         = input_data["initial level of refinement"];
    const double & REDCDT                                   = input_data["CFL condition"];
    const double & T                                        = input_data["final time in seconds"];
    const double & SPACE_ADAPTDT                            = input_data["space adaptation procedure interval in seconds"];
    const double & SAVEDT                                   = input_data["saving interval in seconds"];
    const double & DELTAT                                   = input_data["maximum time step allowed"];
    const double & initial_delta_t 			    = input_data["initial time step"];
    const double & mesh_size_dry                            = input_data["desired resolution in meters of the mesh size in dry regions"];
    const double & mesh_size_wet                            = input_data["minimum resolution in meters of the mesh size in wet regions"];
    const double & mesh_size_interface                      = input_data["desired resolution in meters of the mesh size in wet-dry interface regions"];
    const bool   & is_time_adaptivity                       = input_data["do you want the time step predictor?"];
    const bool   & is_initial_refinement                    = input_data["do you want to refine the mesh initially?"];
    const bool   & is_space_adaptivity                      = input_data["do you want the space adaptation with interface tracking?"];
    const bool   & is_non_reflBC                            = input_data["do you want non reflecting BC?"];
    const bool   & is_isothermal                            = input_data["do you want an isothermal simulation?"];
    const bool   & is_limiter                               = input_data["do you want the FCT limiting?"];
    const bool   & is_max_time_step_from_CFL                = input_data["do you want the maximum time step given by CFL condition for the transport term?"];
    h_min                                                   = input_data["minimum material height threshold"];
    const double & grav                                     = input_data["gravitational field"];
    const double & density                                  = input_data["material density"];
    const double & T_env                                    = input_data["environment temperature"];
    const double & specific_heat_pressure                   = input_data["specific heat at constant pressure"];
    const double & convective_coeff                         = input_data["convection coefficient"];
    const double & tolerance_space_adapt                    = input_data["tolerance space adaptation"];
    const double & sigma_vent                               = input_data["area discrete vent"];
    const double & mesh_size_vent 			    = input_data["mesh size vent"];
    const double & saturation_coeff                         = input_data["saturation coefficient"];
    const double & x_v                                      = input_data["x vent location"];
    const double & y_v                                      = input_data["y vent location"];
    const double & Q_vent                                   = input_data["lava vent discharge"]; 
    const double & T_vent                                   = input_data["lava vent effusion temperature"];
    const double & b_exp_coeff                              = input_data["b coefficient"];
    const double & T_ref                                    = input_data["T_ref"];
    const double & nu_ref                                   = input_data["nu reference"];
    const std::string & SAVE_DIR                            = input_data["home saving directory, i.e., where we can find the directory results"];
    const std::string & DEM_DIR                             = input_data["dem file, complete path"]; 

    L = res*(Nx-1);
    H = res*(Ny-1);

    // Connectivity of local element
    constexpr p4est_topidx_t simple_conn_num_vertices = 4;
    constexpr p4est_topidx_t simple_conn_num_trees = 1;
    const double simple_conn_p[simple_conn_num_vertices*2] =
    {
        0,  0,
        0,  H,
        L,  0,
        L,  H
    };

    const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
    {  1,    3,    4,    2,    1 };

    // Management of solutions ordering
    ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 0> (gt); };
    ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 1> (gt); };
    ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 2> (gt); };
    ordering ordTh = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 3> (gt); };

    // Initialize MPI
    int rank, size;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    // Generate the mesh in 2d
    tmesh tmsh;
    tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
            simple_conn_t, simple_conn_num_trees);

    int recursive = 1;
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (recursive);

    // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
    tmesh::idx_t gn_nodes    = tmsh.num_global_nodes (); // Return total number of nodes owned by all process
    tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes (); // Return number of nodes owned by local process
    tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();  // Return number of quadrants owned by local process across all trees
    tmesh::idx_t gn_elements = tmsh.num_global_quadrants (); // Return number of quadrants owned by all processes across all trees

    // Allocate initial data container
    Q1 sol  (ln_nodes * 4);
    Q1 incr (ln_nodes * 4);
    sol.get_owned_data  ().assign (sol.get_owned_data  ().size (), 0.0);
    incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);

    Q1 mass (ln_nodes * 4);
    bim2a_mass_vector (tmsh, mass, ordh);
    bim2a_mass_vector (tmsh, mass, ordUx);
    bim2a_mass_vector (tmsh, mass, ordUy);
    bim2a_mass_vector (tmsh, mass, ordTh);
    mass.assemble ();

    Q0 sol_onehalf (ln_elements * 4);
    sol_onehalf.get_owned_data  ().assign (sol_onehalf.get_owned_data  ().size (), 0.0);

    Q0 Z_onehalf (ln_elements);
    Z_onehalf.get_owned_data  ().assign (Z_onehalf.get_owned_data  ().size (), 0.0);

    std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 4);

    Q1 Z (ln_nodes);
    Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);

    std::string str = ""; 
    char filename[255]="", arr[255]="";

    str = std::string(DEM_DIR);
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);

    octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
    octave_value v;

    octave_io_open (filename, m_in, &m_out);
    octave_load (VARNAME_1, v);
    Matrix M = v.matrix_value ();
    dem.resize (M.numel ());
    std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());


    // Initialize 
    TIC ();
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
            quadrant != tmsh.end_quadrant_sweep ();
            ++quadrant)
    {
        double xx_c=quadrant->centroid(0);
        double yy_c=quadrant->centroid(1);

        for (int ii = 0; ii < 4; ++ii)
        {
            if (! quadrant->is_hanging (ii)){
                double xx=quadrant->p(0,ii);
                double yy=quadrant->p(1,ii); 

                sol [ordh     (quadrant->gt (ii))] = h0_fun  (xx, yy, grav);
                sol [ordUx    (quadrant->gt (ii))] = Ux0_fun (xx, yy, grav);
                sol [ordUy    (quadrant->gt (ii))] = Uy0_fun (xx, yy, grav);
                sol [ordTh    (quadrant->gt (ii))] = Th0_fun (xx, yy);

                Z   [quadrant->gt (ii)] = dem_fun(xx,yy); 
            }

            else
            {
                // touch parent nodes to set up distributed vector structure
                sol [ordh   (quadrant->gparent(0,ii))] += 0.;
                sol [ordh   (quadrant->gparent(1,ii))] += 0.;
                sol [ordUx  (quadrant->gparent(0,ii))] += 0.;
                sol [ordUx  (quadrant->gparent(1,ii))] += 0.;
                sol [ordUy  (quadrant->gparent(0,ii))] += 0.;
                sol [ordUy  (quadrant->gparent(1,ii))] += 0.;
                sol [ordTh  (quadrant->gparent(0,ii))] += 0.;
                sol [ordTh  (quadrant->gparent(1,ii))] += 0.;

                Z   [quadrant->gparent(0,ii)] += 0.;
                Z   [quadrant->gparent(1,ii)] += 0.;
            }
        }
    }

    // bim2a_solution_with_ghosts in quad_operators.cpp
    bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy, false);
    bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordTh);

    bim2a_solution_with_ghosts (tmsh, Z, replace_op);

    bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy, false);
    bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordTh);


    if (is_initial_refinement)
    {

        Q1 only_h (ln_nodes);
        bim2a_solution_with_ghosts (tmsh, only_h);
        for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
        {
            only_h(idx) = sol(ordh(idx));
        }
        only_h.assemble (replace_op);

        auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);

        auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
        {
            return estimator_grad(q, dh, only_h);
        };
        auto estimator_flux = [& only_h] (tmesh::quadrant_iterator q)
        {
            std::array<double,4> h_mesh = {0,0,0,0};
            for (int ii = 0; ii < 4; ++ii)
            {
                if (! q->is_hanging (ii)){
                    h_mesh[ii] = only_h[q->gt (ii)];
                }
                else
                {
                    h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
                }
            }
            const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 
            return (basin_check); 
        };

        auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
        {
            std::array<double,4> h_mesh = {0,0,0,0};
            for (int ii = 0; ii < 4; ++ii)
            {
                if (! q->is_hanging (ii)){
                    h_mesh[ii] = only_h[q->gt (ii)];
                }
                else
                {
                    h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
                }
            }
            const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])==0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 
            return (basin_check); 
        };

        tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, tolerance_space_adapt, x_v, y_v, mesh_size_vent, 6, 0, 0);
        tmsh.metrics_refine (1e6);  // RAFFINAMENTO (arg is max element)

        // Ottengo i parametri della mesh corrente
        gn_nodes    = tmsh.num_global_nodes ();
        ln_nodes    = tmsh.num_owned_nodes ();
        ln_elements = tmsh.num_local_quadrants ();
        gn_elements = tmsh.num_global_quadrants ();

        // Interpolo sol sulla nuova mesh
        Q1 incr_ (ln_nodes * 4);
        incr_.get_owned_data ().assign (incr_.get_owned_data ().size(), 0.0);
        incr_.assemble ();

        Q1 mass_ (ln_nodes * 4);
        bim2a_mass_vector (tmsh, mass_, ordh );
        bim2a_mass_vector (tmsh, mass_, ordUx);
        bim2a_mass_vector (tmsh, mass_, ordUy);
        bim2a_mass_vector (tmsh, mass_, ordTh);
        mass_.assemble ();

        Q0 sol_onehalf_ (ln_elements * 4);
        sol_onehalf_.get_owned_data ().assign (sol_onehalf_.get_owned_data ().size(), 0.0);

        Q0 Z_onehalf_ (ln_elements);
        Z_onehalf_.get_owned_data ().assign (Z_onehalf_.get_owned_data ().size(), 0.0);

        std::vector<std::array<double,4>> incr_anti_diff_ (ln_elements * 4);

        Q1 sol_ (ln_nodes * 4);
        Q1 Z_ (ln_nodes);
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep ();
                ++quadrant)
        {
            double xx_c=quadrant->centroid(0);
            double yy_c=quadrant->centroid(1);

            for (int ii = 0; ii < 4; ++ii)
            {
                if (! quadrant->is_hanging (ii)){
                    double xx=quadrant->p(0,ii);
                    double yy=quadrant->p(1,ii);

                    sol_ [ordh     (quadrant->gt (ii))] = h0_fun  (xx, yy, grav);
                    sol_ [ordUx    (quadrant->gt (ii))] = Ux0_fun (xx, yy, grav);
                    sol_ [ordUy    (quadrant->gt (ii))] = Uy0_fun (xx, yy, grav);
                    sol_ [ordTh    (quadrant->gt (ii))] = Th0_fun (xx, yy);

                    Z_[quadrant->gt (ii)] = dem_fun(xx,yy);
                }
                else
                {
                    sol_ [ordh   (quadrant->gparent(0,ii))] += 0.;
                    sol_ [ordh   (quadrant->gparent(1,ii))] += 0.;
                    sol_ [ordUx  (quadrant->gparent(0,ii))] += 0.;
                    sol_ [ordUx  (quadrant->gparent(1,ii))] += 0.;
                    sol_ [ordUy  (quadrant->gparent(0,ii))] += 0.;
                    sol_ [ordUy  (quadrant->gparent(1,ii))] += 0.;
                    sol_ [ordTh  (quadrant->gparent(0,ii))] += 0.;
                    sol_ [ordTh  (quadrant->gparent(1,ii))] += 0.;

                    Z_[quadrant->gparent(0,ii)] += 0.;
                    Z_[quadrant->gparent(1,ii)] += 0.;
                }
            }
        }

        bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordh,  false);
        bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUx, false);
        bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUy, false);
        bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordTh);

        bim2a_solution_with_ghosts (tmsh, Z_, replace_op);

        bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordh,  false);
        bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUx, false);
        bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUy, false);
        bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordTh);


        sol                 = sol_;
        incr                = incr_;
        incr_anti_diff      = incr_anti_diff_;
        mass                = mass_;
        sol_onehalf         = sol_onehalf_;
        Z                   = Z_;
        Z_onehalf           = Z_onehalf_;
    }

    Q1 sol_dyn                 = sol;
    Q1 sold_dyn                = sol;
    Q1 soldd_dyn               = sol;
    Q1 sol_2_dyn               = sol;
    Q1 incr_dyn                = incr;
    Q1 incr_vent_dyn           = incr;
    Q1 incr_second_dyn         = incr;
    Q1 P_plus_dyn              = incr;
    Q1 P_minus_dyn             = incr;  
    Q1 mass_dyn                = mass;
    Q1 Z_dyn                   = Z;
    Q0 sol_onehalf_dyn         = sol_onehalf;
    Q0 Z_onehalf_dyn           = Z_onehalf;

    std::vector<std::array<double,4>> incr_anti_diff_dyn = incr_anti_diff;

    TG2_scheme stp(sol_dyn, 
            sold_dyn, 
            soldd_dyn, 
            sol_2_dyn,
            incr_dyn,
            incr_vent_dyn,
            incr_second_dyn,
            incr_anti_diff_dyn,
            P_plus_dyn, 
            P_minus_dyn, 
            sol_onehalf_dyn, 
            mass_dyn,
            ordh, 
            ordUx, 
            ordUy, 
            ordTh,
            Z_dyn,
            Z_onehalf_dyn,
            DELTAT, 
            h_min, 
            is_non_reflBC, 
            is_isothermal, 
            is_limiter,
            grav, 
            nu_ref, 
            T_ref, 
            b_exp_coeff,
            saturation_coeff,
            density, 
            T_env,
            specific_heat_pressure,
            convective_coeff,
            x_v, 
            y_v, 
            Q_vent, 
            T_vent, 
            sigma_vent);


    // Create the results folder if it doesn't exist
    std::string results_dir = std::string(SAVE_DIR) + "/results";
    std::filesystem::create_directories(results_dir);

    // Save initial conditions
    str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d"; 
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);
    tmsh.octbin_export (filename, sol_dyn, ordh);

    str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);
    tmsh.octbin_export (filename, sol_dyn, ordUx); 

    str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);
    tmsh.octbin_export (filename, sol_dyn, ordUy);

    str = std::string(SAVE_DIR) + "/results/swe_Th_%4.4d";
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);
    tmsh.octbin_export (filename, sol_dyn, ordTh);

    str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0); 
    tmsh.octbin_export (filename, Z_dyn);

    std::vector<double> full_time_vector;
    full_time_vector.reserve (static_cast<int> (T/DELTAT));
    std::vector<double> save_time_vector;
    save_time_vector.reserve (static_cast<int> (T/SAVEDT));

    // Time loop
    double time      = 0.0;
    double time_old  = 0.0;
    double time_oldd = 0.0;

    stp.set_dt (DELTAT);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
            quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {

        stp.compute_dt(quadrant);

    }
    double max_dt = REDCDT * stp.dt; 

    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&max_dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);
    stp.set_dt(max_dt);
    stp.set_old_dt(0.);
    time_old  -= stp.dt;
    time_oldd -= 2*stp.dt;
    stp.set_times(time, time_old, time_oldd);  

    if (rank==0) {

        full_time_vector.push_back (0.0);
        save_time_vector.push_back (0.0);

    }
    int count = 0;

    double savecount = 0.0, space_adapt_count = 0.0;

    if (rank == 0)
    {
        std::cout << "start loop" << std::endl;
    }

    int counter_savings = 0, tot_number_savings = std::round(T/SAVEDT);

    TIC();
    while (counter_savings != tot_number_savings)
    {
        // Reset increment, and limiter terms
        incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
        incr_dyn.assemble (replace_op);

        incr_second_dyn.get_owned_data ().assign (incr_second_dyn.get_owned_data ().size (), 0.0);
        incr_second_dyn.assemble (replace_op);

        P_plus_dyn.get_owned_data ().assign (P_plus_dyn.get_owned_data ().size (), 0.0);
        P_plus_dyn.assemble (replace_op);

        P_minus_dyn.get_owned_data ().assign (P_minus_dyn.get_owned_data ().size (), 0.0);
        P_minus_dyn.assemble (replace_op);

        // compute time step, 
        stp.Fr = 0.;
        stp.set_dt (DELTAT);
        if (is_max_time_step_from_CFL)
        {
            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                    quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
            {
                stp.compute_dt(quadrant);
            }
        }
        max_dt = REDCDT * stp.dt;

        stp.g_coeff = 1./(1.-stp.Fr*stp.Fr);

        stp.set_dt(time==0 ? std::min(initial_delta_t, max_dt) : max_dt); // deltat max
        MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);

        // Print current time
        if (rank==0)
        {
            std::cout << "MAXIMUM TIME STEP = " << stp.dt << std::endl;
        }   

        // time adaptivity
        if (is_time_adaptivity)
        { 
            stp.nu_htot = 0.;
            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                    quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
            {
                stp.compute_dt_adaptive(quadrant);
            }
            MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.nu_htot), 1, MPI_DOUBLE, MPI_SUM, tmsh.comm);

            static constexpr double local_estimator_time_tolerance = 1e-5;

            const double candidate_dt = local_estimator_time_tolerance/std::sqrt(stp.nu_htot)*(stp.time-stp.timed);
            stp.set_dt( (stp.nu_htot>0 && candidate_dt<stp.dt) ? candidate_dt : stp.dt );
        }

        if (stp.dt == 0 && rank == 0)
        {
            std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
            exit( -1. );
        }

        // check save with given frequency
        stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? SAVEDT-savecount : stp.dt);

        time_oldd = time_old;
        time_old = time;
        time += stp.dt; 
        savecount += stp.dt;
        space_adapt_count += stp.dt;

        // Print current time
        if (rank==0) 
        {
            std::cout << "TIME = " << time << ", dt = " << stp.dt << std::endl;
            full_time_vector.push_back (time);
        }

        // first step computes the q^(n,2)
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
            stp.first_step(quadrant);
        }

        // 
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
            stp.compute_nodal_anti_diffusive_fluxes(quadrant);
        }
        incr_dyn.assemble ();
        incr_second_dyn.assemble ();
        P_plus_dyn.assemble ();
        P_minus_dyn.assemble ();

        stp.set_times(time, time_old, time_oldd);
        soldd_dyn = sold_dyn;
        sold_dyn  = sol_dyn;   

        // low order solution
        for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
        {
            sol_dyn.get_owned_data ()[kk] += stp.dt_expl_32*incr_dyn.get_owned_data ()[kk]/mass_dyn.get_owned_data ()[kk];
        }
        sol_dyn.assemble(replace_op); 

        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep ();
                ++quadrant)
        {
            for (int ii = 0; ii < 4; ++ii)
            {
                if (! quadrant->is_hanging (ii) && sol_dyn [ordh    (quadrant->gt (ii))] < 0) {
                    sol_dyn [ordh     (quadrant->gt (ii))] = 0.; 
                    sol_dyn [ordTh    (quadrant->gt (ii))] = 0.;
                }
            }
        }
        bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordh,  false);
        bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUx, false);
        bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUy, false);
        bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordTh);

        incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
        incr_dyn.assemble (replace_op);
 
        incr_vent_dyn.get_owned_data ().assign (incr_vent_dyn.get_owned_data ().size (), 0.0);
        incr_vent_dyn.assemble (replace_op);

        // second order correction
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
            stp.second_step(quadrant);
        }
        incr_dyn.assemble ();
        incr_vent_dyn.assemble ();

        // Compute here the second step solution!, i.e., q^(n,3)
        for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=4)
        {
            stp.solve_non_lin(kk);
        }
        sol_2_dyn.assemble(replace_op);
        sol_dyn  .assemble(replace_op);

        // solution update
        incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
        incr_dyn.assemble (replace_op);
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
            stp.compute_updated_sol(quadrant);
        }
        incr_dyn.assemble ();

        for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=4)
        {
            stp.compute_updated_sol(kk);
        }
        sol_dyn.assemble (replace_op);

        // Save solution
        if ((savecount-SAVEDT) >= -std::numeric_limits<double>::epsilon()*SAVEDT) 
        {
            if (rank == 0)
                std::cout << "savecount = " << savecount << std::endl;
            count++;
            save_time_vector.push_back (time);

            str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d";
            strcpy(arr, str.c_str());
            sprintf(filename, arr,   count);
            tmsh.octbin_export (filename, sol_dyn, ordh);

            str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
            strcpy(arr, str.c_str());
            sprintf(filename, arr,  count);
            tmsh.octbin_export (filename, sol_dyn, ordUx);

            str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
            strcpy(arr, str.c_str());
            sprintf(filename, arr,  count);
            tmsh.octbin_export (filename, sol_dyn, ordUy);

            str = std::string(SAVE_DIR) + "/results/swe_Th_%4.4d";
            strcpy(arr, str.c_str());
            sprintf(filename, arr,  count);
            tmsh.octbin_export (filename, sol_dyn, ordTh);

            str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
            strcpy(arr, str.c_str());
            sprintf(filename, arr,  count);
            tmsh.octbin_export (filename, Z_dyn);
            savecount = 0.0;

            counter_savings++;
        }


        if (is_space_adaptivity &&  ((space_adapt_count-SPACE_ADAPTDT) >= -std::numeric_limits<double>::epsilon()*SPACE_ADAPTDT))
        {
            Q1 only_h (ln_nodes);
            bim2a_solution_with_ghosts (tmsh, only_h);
            for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
            {
                only_h(idx) = sol_dyn(ordh(idx));
            }
            only_h.assemble (replace_op);

            Q1 only_Ux (ln_nodes);
            bim2a_solution_with_ghosts (tmsh, only_Ux);
            for (auto idx = only_Ux.get_range_start (); idx != only_Ux.get_range_end (); ++idx)
            {
                only_Ux(idx) = sol_dyn(ordUx(idx));
            }
            only_Ux.assemble (replace_op);

            Q1 only_Uy (ln_nodes);
            bim2a_solution_with_ghosts (tmsh, only_Uy);
            for (auto idx = only_Uy.get_range_start (); idx != only_Uy.get_range_end (); ++idx)
            {
                only_Uy(idx) = sol_dyn(ordUy(idx));
            }
            only_Uy.assemble (replace_op);
       
            std::set<int> global_index_quad;
            for (auto q = tmsh.begin_quadrant_sweep ();
                    q != tmsh.end_quadrant_sweep ();
                    ++q)
            {
                quadrant_marker_list(q, only_h,  only_Ux, only_Uy, stp.dt, global_index_quad);
            }      


            auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);
            auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
            {
                return estimator_grad(q, dh, only_h);
            };
  
            auto estimator_flux = [& global_index_quad] (tmesh::quadrant_iterator q)
            {
                if ( global_index_quad.find(q->get_global_quad_idx ()) != global_index_quad.end() )
                {
                    return 1;
                }
                return 0;
            };

            auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
            {
                std::array<double,4> h_mesh = {0,0,0,0};
                for (int ii = 0; ii < 4; ++ii)
                {
                    if (! q->is_hanging (ii)){
                        h_mesh[ii] = only_h[q->gt (ii)];
                    }
                    else
                    {
                        h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
                    }
                }
                const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])<h_min && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])<h_min) ? 1 : 0;
                return (basin_check); 
            };

            tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, tolerance_space_adapt, 6, 0, 0);
            tmsh.metrics_refine (1e6);  // RAFFINAMENTO (arg is max element)

            // Ottengo i parametri della mesh corrente
            gn_nodes    = tmsh.num_global_nodes ();
            ln_nodes    = tmsh.num_owned_nodes ();
            ln_elements = tmsh.num_local_quadrants ();
            gn_elements = tmsh.num_global_quadrants ();

            // Interpolate on the new mesh
            Q1 sol (ln_nodes * 4);
            bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
            bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
            bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy, false);
            bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordTh);
            interpolate_vector (tmsh, sol_dyn, sol, ordh);
            interpolate_vector (tmsh, sol_dyn, sol, ordUx);
            interpolate_vector (tmsh, sol_dyn, sol, ordUy);
            interpolate_vector (tmsh, sol_dyn, sol, ordTh);
            sol.assemble (replace_op);

            Q1 sold (ln_nodes * 4);
            bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordh,  false);
            bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUx, false);
            bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUy, false);
            bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordTh);
            interpolate_vector (tmsh, sold_dyn, sold, ordh);
            interpolate_vector (tmsh, sold_dyn, sold, ordUx);
            interpolate_vector (tmsh, sold_dyn, sold, ordUy);
            interpolate_vector (tmsh, sold_dyn, sold, ordTh);
            sold.assemble (replace_op);

            Q1 soldd (ln_nodes * 4);
            bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordh,  false);
            bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUx, false);
            bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUy, false);
            bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordTh);
            interpolate_vector (tmsh, soldd_dyn, soldd, ordh );
            interpolate_vector (tmsh, soldd_dyn, soldd, ordUx);
            interpolate_vector (tmsh, soldd_dyn, soldd, ordUy);
            interpolate_vector (tmsh, soldd_dyn, soldd, ordTh);
            soldd.assemble (replace_op);

            Q1 incr (ln_nodes * 4);
            incr.get_owned_data ().assign (incr.get_owned_data ().size(), 0.0);
            incr.assemble ();

            std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 4);

            Q1 mass (ln_nodes * 4);
            bim2a_mass_vector (tmsh, mass, ordh );
            bim2a_mass_vector (tmsh, mass, ordUx);
            bim2a_mass_vector (tmsh, mass, ordUy);
            bim2a_mass_vector (tmsh, mass, ordTh);
            mass.assemble ();

            Q0 sol_onehalf (ln_elements * 4);
            sol_onehalf.get_owned_data ().assign (sol_onehalf.get_owned_data ().size(), 0.0);

            Q0 Z_onehalf (ln_elements);
            Z_onehalf.get_owned_data ().assign (Z_onehalf.get_owned_data ().size(), 0.0);

            Q1 Z (ln_nodes);
            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                    quadrant != tmsh.end_quadrant_sweep ();
                    ++quadrant)
            {
                double xx_c=quadrant->centroid(0);
                double yy_c=quadrant->centroid(1); 
                for (int ii = 0; ii < 4; ++ii)
                {
                    if (! quadrant->is_hanging (ii)){
                        double xx=quadrant->p(0,ii);
                        double yy=quadrant->p(1,ii);
                        Z[quadrant->gt (ii)] = dem_fun(xx,yy); 
                    }
                    else
                    {
                        Z[quadrant->gparent(0,ii)] += 0.;
                        Z[quadrant->gparent(1,ii)] += 0.;
                    }
                }
            }
            bim2a_solution_with_ghosts (tmsh, Z,    replace_op);
            bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
            bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
            bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy, false);
            bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordTh);

            sol_dyn                 = sol;
            sold_dyn                = sold;
            soldd_dyn               = soldd;
            sol_2_dyn               = soldd;
            incr_dyn                = incr;
            incr_second_dyn         = incr;
            incr_anti_diff_dyn      = incr_anti_diff;
            P_plus_dyn              = incr;
            P_minus_dyn             = incr;
            mass_dyn                = mass;
            sol_onehalf_dyn         = sol_onehalf;
            Z_onehalf_dyn           = Z_onehalf;
            Z_dyn                   = Z;

            space_adapt_count = 0.0;
        }
    }

    if (rank == 0)
    {
        str = std::string(SAVE_DIR) + "/results/timesteps.octbin"; 
        strcpy(arr, str.c_str());
        sprintf(filename, arr, 0);

        ColumnVector vtmp (save_time_vector.size ());
        std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
        octave_io_mode m;
        octave_io_open (filename, gz_write_mode, &m);
        octave_save ("save_time", vtmp);
        vtmp.resize (full_time_vector.size ());
        std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
        octave_save ("full_time", vtmp);
        octave_io_close ();
    }

    TOC ("loop completed");

    // Close MPI and print report
    MPI_Barrier (MPI_COMM_WORLD);
    if (rank == 0) { print_timing_report (); }
    MPI_Finalize ();
    return 0;
}
