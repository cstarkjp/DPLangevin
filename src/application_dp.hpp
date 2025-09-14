// 
// "Dornic" operator-splitting method of integrating DP-type Langevin equations
//   adapted from code by Paula Villa Martín and Victor Buendía.
//
// CPS 2025-09-02
// 

#ifndef APPLICATION_DP_HPP
#define APPLICATION_DP_HPP

namespace py = pybind11;
typedef py::array_t<double, py::array::c_style> py_array_t;

class Results 
{
private:
    const int n_epochs, n_cells, n_x, n_y, n_z;
    py_array_t return_epochs, return_mean_densities, density;

public:
    Results(
        const int n_epochs, const int n_cells, 
        const int n_x, const int n_y, const int n_z
    ) : n_epochs(n_epochs), n_cells(n_cells), n_x(n_x), n_y(n_y), n_z(n_z) {}

    void prep_epochs(const dbl_vec_t& epochs)
    {
        py_array_t epochs_array(n_epochs);
        auto epochs_proxy = epochs_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            epochs_proxy(i) = epochs[i];
        };
        return_epochs = epochs_array;
    }
    py_array_t get_epochs() const { return return_epochs; }

    void prep_mean_densities(const dbl_vec_t& mean_densities)
    {
        py_array_t mean_densities_array(n_epochs);
        auto mean_densities_proxy = mean_densities_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            mean_densities_proxy(i) = mean_densities[i];
        };
        return_mean_densities = mean_densities_array;
    }
    py_array_t get_mean_densities() const { return return_mean_densities; }

    void prep_density(const dbl_vec_t& cell_density)
    {
        py_array_t density_array({n_x, n_y});
        auto density_proxy = density_array.mutable_unchecked();
        // ASSERT: n_cells = n_x * n_y * n_z
        for (auto i=0; i<n_cells; i++)
        {
            // density_proxy(i, 0) = epochs[i];
            // density_proxy(i, 1) = mean_densities[i];
        };
        density = density_array;
    }
    py_array_t get_density() const { return density; }
};

auto dp(
    const double linear, const double quadratic, 
    const double diffusion, const double noise, 
    const double t_max, const double dx, const double dt, 
    const int random_seed,
    const GridDimension grid_dimension, 
    const int_vec_t& grid_size,
    const GridTopology grid_topology,
    const BoundaryCondition boundary_condition,
    const InitialCondition initial_condition,
    const IntegrationMethod integration_method
) -> Results;

class DPLangevin : public LangevinBase 
{
public:
    double quadratic_coeff;
    double D;

    DPLangevin() = default;
    DPLangevin(Parameters params)  
    {
        dt = params.dt;
        dx = params.dx;
        dtm = 0.5*dt;
        dts = dt/6.0;
        n_cells = params.n_cells;
        cell_density = dbl_vec_t(n_cells, 0.0); 
        aux_cell_new = dbl_vec_t(n_cells);
        aux_cell_old = dbl_vec_t(n_cells);
        k1 = dbl_vec_t(n_cells, 0.0);
        k2 = dbl_vec_t(n_cells, 0.0);
        k3 = dbl_vec_t(n_cells, 0.0);
        k4 = dbl_vec_t(n_cells, 0.0);
    }

    void set_nonlinear_coefficients(const Coefficients &f_coefficients)
    override
    {
        quadratic_coeff = f_coefficients.quadratic;
        D = f_coefficients.diffusion / (dx*dx);
    }

    void check()
    {
        std::cout << "DPLangevin::  check"  << std::endl;
    }

    auto nonlinear_rhs(const int i_cell, const dbl_vec_t &field) 
    const -> double
    override
    {
        // Non-linear terms
        const double quadratic_term 
            = -quadratic_coeff*field[i_cell]*field[i_cell];

        // Integration of diffusion
        double diffusion_sum = 0.0;
        int n_neighbors = neighbors[i_cell].size();
        for (auto i=0; i<n_neighbors; i++)
        {
            auto i_neighbor = neighbors[i_cell][i];
            diffusion_sum += field[i_neighbor];
        }
        diffusion_sum = D*(diffusion_sum - n_neighbors*field[i_cell]);
        return diffusion_sum + quadratic_term;
    }
};

class SimDP 
{
private:
    Coefficients f_coeffs;
    Parameters p;
    RNG rng; 
    DPLangevin *dpLangevin;
    int n_epochs;
    py_array_t return_epochs, return_mean_densities, density;
    dbl_vec_t epochs;
    dbl_vec_t mean_densities;
    bool is_initialized = false;

    void construct_grid()
    {
        std::cout << "construct_grid::  dpLangevin = " << dpLangevin << std::endl;
        switch (p.grid_dimension)
        {
            case (GridDimension::D1):
                dpLangevin->construct_1D_grid(p);
                break;
            case (GridDimension::D2):
            default:
                dpLangevin->construct_2D_grid(p);
                break;
        }    
    }

    void initialize_grid()
    {
        switch (p.initial_condition)
        {
            case (InitialCondition::RANDOM_GAUSSIAN):
                dpLangevin->ic_random_uniform(rng);
                break;
            case (InitialCondition::CONSTANT_VALUE):
                dpLangevin->ic_constant_value(1.0);
                break;
            case (InitialCondition::SINGLE_SEED):
                dpLangevin->ic_single_seed(p.n_cells/2, 1.0);
                break;
            case (InitialCondition::RANDOM_UNIFORM):
            default:
                dpLangevin->ic_random_uniform(rng);
                break;
        }  
    }

    int count_epochs()
    {
        // Count total number of time steps, 
        //    just in case rounding causes problems
        int n_epochs;
        double t; 
        for (
            n_epochs=0, t=0; 
            t<=p.t_max+p.dt; 
            t+=p.dt, n_epochs++
        ) {}
        return n_epochs;
    }

    bool integrate(dbl_vec_t& epochs, dbl_vec_t& mean_densities)
    {
        int i;
        double t; 

        std::cout << "integrate::  n_epochs = " << n_epochs << std::endl;
        std::cout << "integrate::  epochs.size = " << epochs.size() << std::endl;
        std::cout << "integrate::  mean_densities.size = " << mean_densities.size() << std::endl;

        switch (p.integration_method)
        {
            case (IntegrationMethod::EULER):
                std::cout << "integrate::  Euler "<< std::endl;
                for (i=0, t=0; i<epochs.size(); t+=p.dt, i++)
                {
                    dpLangevin->integrate_euler(rng);
                    epochs[i] = t;
                    mean_densities[i] = dpLangevin->get_mean_density();
                };
                break;
            case (IntegrationMethod::RUNGE_KUTTA):
            default:
                std::cout << "integrate::  Runge-Kutta "<< std::endl;
                // std::cout << "integrate::  rng" << rng << std::endl;
                std::cout << "integrate::  dpLangevin = " << dpLangevin << std::endl;
                for (i=0, t=0; i<epochs.size(); t+=p.dt, i++)
                {
                    dpLangevin->integrate_rungekutta(rng);
                    epochs[i] = t;
                    mean_densities[i] = dpLangevin->get_mean_density();
                };
                break;
        }
        return true;
    }

    bool prep_epochs()
    {
        py_array_t epochs_array(n_epochs);
        auto epochs_proxy = epochs_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            epochs_proxy(i) = epochs[i];
        };
        return_epochs = epochs_array;
        return true;
    }

    bool prep_mean_densities()
    {
        py_array_t mean_densities_array(n_epochs);
        auto mean_densities_proxy = mean_densities_array.mutable_unchecked();
        for (auto i=0; i<n_epochs; i++)
        {
            mean_densities_proxy(i) = mean_densities[i];
        };
        return_mean_densities = mean_densities_array;
        return true;
    }

    bool prep_density()
    {
        py_array_t density_array({p.n_x, p.n_y});
        auto density_proxy = density_array.mutable_unchecked();
        // std::cout << "prep_density: p.n_cells = " << p.n_cells << std::endl;
        // std::cout << "prep_density: p.n_x = " << p.n_x << std::endl;
        // std::cout << "prep_density: p.n_y = " << p.n_y << std::endl;
        // std::cout << "prep_density: p.n_z = " << p.n_z << std::endl;
        if (not (p.n_cells == p.n_x * p.n_y * p.n_z)) { 
            std::cout << "prep_density: failed" << std::endl;
            return false; 
        }
        for (auto i=0; i<p.n_cells; i++)
        {
            // density_proxy(i, 0) = epochs[i];
            // density_proxy(i, 1) = mean_densities[i];
        };
        density = density_array;
        return true;
    }


public:
    // SimDP() = default;
    SimDP(
        const double linear, const double quadratic,
        const double diffusion, const double noise, 
        const double t_max, const double dx, const double dt, 
        const int random_seed,
        const GridDimension grid_dimension,
        const int_vec_t& grid_size,
        const GridTopology grid_topology,
        const BoundaryCondition boundary_condition,
        const InitialCondition initial_condition,
        const IntegrationMethod integration_method
    ) : f_coeffs(linear, quadratic, diffusion, noise),
        p(
            t_max, dx, dt, random_seed,
            grid_dimension, grid_size, grid_topology, 
            boundary_condition, initial_condition, integration_method
        )
    {
        RNG rng(p.random_seed); 
        dpLangevin = new DPLangevin(p);
        f_coeffs.print();
        p.print();
    }

    bool initialize()
    {
        construct_grid();
        initialize_grid();
        dpLangevin->check();
        dpLangevin->set_coefficients(f_coeffs);
        is_initialized = true;
        return true;
    }

    bool run(void)
    {
        if (not is_initialized) 
        { 
            std::cout << "Failure: must initialize first" << std::endl;
            return false; 
        }
        n_epochs = count_epochs();
        epochs = dbl_vec_t(n_epochs, 0.0);
        mean_densities = dbl_vec_t(n_epochs, 0.0);
        bool did_integrate = integrate(epochs, mean_densities);
        bool did_finalize = (
            did_integrate and
            prep_epochs() and prep_mean_densities() and prep_density()
        ); 
        return did_finalize;
    }

    py_array_t get_epochs() const { return return_epochs; }
    py_array_t get_mean_densities() const { return return_mean_densities; }
    py_array_t get_density() const { return density; }
};


#endif
