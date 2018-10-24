#ifndef ADUPROP_MC_HPP_

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include "alg.hpp"
#include "ad.hpp"

extern System sys;

class mc {
    private:
        const gsl_rng_type *T;
        gsl_rng *r;
        size_t mc_samples;
        gsl_vector *gsl_mu;
        gsl_matrix *gsl_cov;
        gsl_matrix *gsl_samples_mat;
        gsl_vector **gsl_samples;
        pVector<double> *samples;
        size_t dim;
        System &sys;
    public:
    mc(size_t mc_samples_, System &sys_) : mc_samples(mc_samples_), sys(sys_) {
        dim = sys.dim();
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        gsl_mu = gsl_vector_alloc (dim);
        gsl_cov = gsl_matrix_alloc(dim, dim);
        gsl_samples_mat = gsl_matrix_alloc(mc_samples, dim);
        gsl_samples = new gsl_vector*[mc_samples];
        samples = new pVector<double>[mc_samples];
        for (size_t i = 0; i < mc_samples; i++)
        {
            gsl_samples[i] = gsl_vector_alloc(dim);
            samples[i].alloc(dim);
        }
    }
    void integrate(pVector<double> &mu0, pMatrix<double> &cv0, size_t tsteps) {
        /* Random number generator */
        std::cout << "Starting MC" << std::endl;
        pVector<double> mu(dim);
        pMatrix<double> cov(dim, dim);
        ad drivers(sys);
        for (size_t i = 0; i < dim; i++)
        {
            gsl_vector_set(gsl_mu, i, mu0[i]);
            for (size_t j = 0; j < dim; j++)
            {
                gsl_matrix_set(gsl_cov, i, j, cv0[i][j]);
            }
        }
        gsl_linalg_cholesky_decomp(gsl_cov);
        for (size_t i = 0; i < mc_samples; i++)
        {
            if (gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_cov, gsl_samples[i]))
            {
                std::cout << "Gaussian input sampling failed" << std::endl;
            }
            else
            {
                for (size_t j = 0; j < dim; j++)
                {
                    samples[i][j] = gsl_vector_get(gsl_samples[i], j);
                }
            }
        }
        for (size_t t = 0; t < tsteps; ++t)
        {
            for (size_t i = 0; i < mc_samples; i++)
            {
                drivers.integrate(samples[i]);
                for (size_t j = 0; j < dim; j++)
                {
                    gsl_vector_set(gsl_samples[i], j, samples[i][j]);
                }
            }
            for (size_t i = 0; i < mc_samples; i++)
            {
                for (size_t j = 0; j < dim; j++)
                {
                    gsl_matrix_set(gsl_samples_mat, i, j, gsl_vector_get(gsl_samples[i], j));
                }
            }
            if (gsl_ran_multivariate_gaussian_mean(gsl_samples_mat, gsl_mu))
            {
                std::cout << "Error in computing mean." << std::endl;
                return;
            }
            if (gsl_ran_multivariate_gaussian_vcov(gsl_samples_mat, gsl_cov))
            {
                std::cout << "Error in computing covariance." << std::endl;
                return;
            }
            for (size_t i = 0; i < dim; i++)
            {
                mu[i] = gsl_vector_get(gsl_mu, i);
                for (size_t j = 0; j < dim; ++j)
                {
                    cov[i][j] = gsl_matrix_get(gsl_cov, i, j);
                }
            }
    }
    for (size_t i = 0; i < dim; i++)
    {
      std::cout << mu[i] << " ";
    }
    cout << std::endl;
    for (size_t i = 0; i < dim; i++)
    {
      std::cout << cov[i][i] << " ";
    }
    cout << std::endl;

    }
    ~mc(){
        gsl_rng_free(r);
        gsl_vector_free(gsl_mu);
        gsl_matrix_free(gsl_cov);
        gsl_matrix_free(gsl_samples_mat);
        delete[] samples;
        for (size_t i = 0; i < mc_samples; i++)
        {
        gsl_vector_free(gsl_samples[i]);
        }
    }

};

#endif  // ADUPROP_AD_HPP_