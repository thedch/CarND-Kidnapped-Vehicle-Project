/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *  Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>

#include "particle_filter.h"

#define YAW_RATE_THRESHOLD 0.0001

using namespace std;

std::vector<LandmarkObs> transform_observations(Particle current_location,
    const std::vector<LandmarkObs> &observations);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first
    // position (based on estimates of x, y, theta and their uncertainties
    // from GPS) and all weights to 1. Add random Gaussian noise to each
    // particle. NOTE: Consult particle_filter.h for more information about this
    // method (and others in this file)

    normal_distribution<double> x_dist{x, std[0]}; // mean, std_dev
    normal_distribution<double> y_dist{y, std[1]};
    normal_distribution<double> theta_dist{theta, std[2]};

    num_particles = 100;

    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.id = i;
        p.x = x_dist(gen);
        p.y = y_dist(gen);
        p.theta = theta_dist(gen);
        p.weight = 1;

        particles.push_back(p);
        weights.push_back(p.weight);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity,
    double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and
    // std::default_random_engine useful.
    // http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    // http://www.cplusplus.com/reference/random/default_random_engine

    double new_x = 0;
    double new_y = 0;
    double new_theta = 0;

    for (auto &particle : particles) {
        // Handle the case where yaw rate is very small
        if (abs(yaw_rate) < YAW_RATE_THRESHOLD) {
            new_x = particle.x + velocity * delta_t * cos(particle.theta);
            new_y = particle.y + velocity * delta_t * sin(particle.theta);
            new_theta = particle.theta;
        } else { // General case where yaw rate is sufficiently large
            new_x = particle.x + velocity/yaw_rate * ( sin(particle.theta +
                (yaw_rate*delta_t)) - sin(particle.theta));

            new_y = particle.y + velocity/yaw_rate * (-cos(particle.theta +
                (yaw_rate*delta_t)) + cos(particle.theta));

            new_theta = particle.theta + yaw_rate * delta_t;
        }

        normal_distribution<double> x_dist{new_x, std_pos[0]}; // mean, std_dev
        normal_distribution<double> y_dist{new_y, std_pos[1]};
        normal_distribution<double> theta_dist{new_theta, std_pos[2]};

        particle.x = x_dist(gen);
        particle.y = y_dist(gen);
        particle.theta = theta_dist(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<single_landmark_s> nearby_landmarks,
    std::vector<LandmarkObs>& transformed_obs, double max_dist) {
    // Puts the index of each nearby_landmarks nearest to each transformed_obs
    // in the ID field of the transformed_obs element. This tells us how far off
    // each measurement is, which is needed for particle weighting
    for (LandmarkObs &obs : transformed_obs) {
        obs.dist_to_landmark = max_dist + 1;
        for (single_landmark_s landmark : nearby_landmarks) {
            double new_dist = dist(obs.x, obs.y, landmark.x_f, landmark.y_f);
            if (new_dist < obs.dist_to_landmark) {
                // Save which landmark is nearest, as well as the distance to that landmark
                obs.id = landmark.id_i;
                obs.dist_to_landmark = new_dist;
                obs.x_dist_to_landmark = abs(obs.x - landmark.x_f);
                obs.y_dist_to_landmark = abs(obs.y - landmark.y_f);
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

    // Set up standard deviation measurements beforehand for simplicity
    float x_std = std_landmark[0];
    float y_std = std_landmark[1];
    float x_std_2 = pow(x_std, 2);
    float y_std_2 = pow(y_std, 2);

    for (Particle &particle : particles) {
        // Build a list of landmarks within range of the particle
        std::vector<single_landmark_s> nearby_landmarks;
        for (single_landmark_s landmark : map_landmarks.landmark_list) {
            if (dist(particle.x, particle.y, landmark.x_f, landmark.y_f) < sensor_range) {
                nearby_landmarks.push_back(landmark);
            }
        }

        // Convert all observations to global frame based on current particle location
        std::vector<LandmarkObs> transformed_obs = transform_observations(particle, observations);
        dataAssociation(nearby_landmarks, transformed_obs, sensor_range);

        // For each observation, use the measurement error to calculate particle
        // weight using multivariate gaussian equation
        float new_weight = 1;
        for (auto obs : transformed_obs) {
            new_weight *= 1/(2*M_PI*x_std*y_std) *
                exp(-(pow(obs.x_dist_to_landmark, 2)/(2*x_std_2) + pow(obs.y_dist_to_landmark, 2)/(2*y_std_2)));
        }
        particle.weight = new_weight;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional
    // to their weight. NOTE: You may find std::discrete_distribution helpful
    // here.
    // http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // Update the weights vector as a convenience
    for (int i = 0; i < num_particles; i++) {
        weights[i] = particles[i].weight;
    }

    discrete_distribution<int> distribution(weights.begin(), weights.end());

    std::vector<Particle> resample_particles;

    for (int i = 0; i < num_particles; i++ ) {
        resample_particles.push_back(particles[distribution(gen)]);
    }

    particles = resample_particles;
}

void ParticleFilter::SetAssociations(Particle& particle,
    const std::vector<int>& associations,
    const std::vector<double>& sense_x,
    const std::vector<double>& sense_y) {
    // particle: the particle to assign each listed association, and
    // association's (x,y) world coordinates mapping to associations: The
    // landmark id that goes along with each listed association sense_x: the
    // associations x mapping already converted to world coordinates sense_y:
    // the associations y mapping already converted to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1); // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1); // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1); // get rid of the trailing space
    return s;
}

// The first task is to transform each observation marker from the vehicle's
// coordinates to the map's coordinates, with respect to our particle.
std::vector<LandmarkObs> transform_observations(Particle current_location,
    const std::vector<LandmarkObs> &observations) {

    std::vector<LandmarkObs> transformed_obs;

    double x = current_location.x;
    double y = current_location.y;
    double theta = current_location.theta;

    for (auto obs : observations) {
        LandmarkObs tfm_obs;
        tfm_obs.x = x + cos(theta)*obs.x - sin(theta)*obs.y;
        tfm_obs.y = y + sin(theta)*obs.x + cos(theta)*obs.y;
        transformed_obs.push_back(tfm_obs);
    }
    return transformed_obs;
}
