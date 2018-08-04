/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *  Author: Tiffany Huang
 */

#include <random>
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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first
    // position (based on estimates of x, y, theta and their uncertainties
    // from GPS) and all weights to 1. Add random Gaussian noise to each
    // particle. NOTE: Consult particle_filter.h for more information about this
    // method (and others in this file)

    default_random_engine gen;

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

    default_random_engine gen;

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

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
    std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed
    // measurement and assign the observed measurement to this particular
    // landmark. NOTE: this method will NOT be called by the grading code. But
    // you will probably find it useful to implement this method and use it as
    // a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian
    // distribution. You can read more about this distribution here:
    // https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your
    // particles are located according to the MAP'S coordinate system. You will
    // need to transform between the two systems. Keep in mind that this transformation
    // requires both rotation AND translation (but no scaling). The following is
    // a good resource for the theory:
    // https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    // and the following is a good resource for the actual equation to implement
    // (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html


    // 1: convert obs to map reference frame
    Particle current_location = get_best_particle(particles);
    transform_observations(current_location, observations);

    // 2: Compare all observations to all landmarks
    for (auto obs : observations) {
        float smallest_dist = std::numeric_limits<int>::max();
        for (auto landmark : map_landmarks.landmark_list) {
            if (dist(obs.x, obs.y, landmark.x_f, landmark.x_f) < smallest_dist) {
                smallest_dist = dist(obs.x, obs.y, landmark.x_f, landmark.x_f);
                obs.id = landmark.id_i;
            }
        }
    }

    for (auto particle : particles) {
        // update weights ???

    }



    // Iterate through all landmark observations, compare with map landmarks
    // There are 42 map landmarks
    // Use nearest neighbor to figure out which obs is which landmark
    // Calculate location as the avg of the obs?

}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional
    // to their weight. NOTE: You may find std::discrete_distribution helpful
    // here.
    // http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;
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

// Gets the best particle, which will be used as the car's current location
Particle get_best_particle(const std::vector<Particle>& particles) {
    int best_weight = 0;
    Particle best;
    for (auto &p : particles) {
        if (p.weight > best_weight) {
            best = p;
            best_weight = best.weight;
        }
    }
    return best;
}

// The first task is to transform each observation marker from the vehicle's
// coordinates to the map's coordinates, with respect to our particle.
void transform_observations(Particle current_location,
    std::vector<LandmarkObs> &observations) {

    double x = current_location.x;
    double y = current_location.y;
    double theta = current_location.theta;

    for (auto obs : observations) {
        obs.x = x + cos(theta)*obs.x - sin(theta)*obs.y;
        obs.y = y + sin(theta)*obs.x + cos(theta)*obs.y;
    }
}
