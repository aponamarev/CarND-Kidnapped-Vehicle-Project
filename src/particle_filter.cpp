/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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

#include "particle_filter.h"
#define EPS 0.0001f

using namespace std;

// Create random number engine
const default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // Set the number of particles.
    num_particles = 1000;
    // Add random Gaussian noise to each particle.
    normal_distribution<double> x_dist(0, std[0]);
    normal_distribution<double> y_dist(0, std[1]);
    normal_distribution<double> theta_dist(0, std[2]);
    // Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    for (int i=0; i<num_particles; i++) {
        Particle particle;
        particle.x = x + x_dist(gen);
        particle.y = y + y_dist(gen);
        particle.theta = theta + theta_dist(gen);
        particle.weight = 1.0f;
        
        particles.push_back(particle);
    }
    
    is_initialized = true;
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // Motion model
    normal_distribution<double> x_dist(0, std_pos[0]);
    normal_distribution<double> y_dist(0, std_pos[1]);
    normal_distribution<double> theta_dist(0, std_pos[2]);
    // Add measurements to each particle
    for (int i=0; i<particles.size(); i++) {
        Particle p = particles[i];
        // f(x) = x(k) + ∫∂x/∂t + U
        // ∫∂x/∂t = ∫​ k​->(k+1​)​​ v * cos(yaw)
        // U = ∫​ k​->(k+1​)​​ v * cos(yaw_dot * ∆t)
        // f(x) = ∫​ k​->(k+1​)​​ v * cos(yaw + yaw_dot * ∆t)
        // f(x) = v/yaw_dot * (sin(yaw+yaw_dot * ∆t) - sin(yaw))
        
        // f(y) = ∫​ k​->(k+1​)​​ v * sin(yaw + yaw_dot * ∆t)
        // f(y) = v/yaw_dot * (cos(yaw) - cos(yaw+yaw_dot * ∆t))
        if (fabs(yaw_rate)>EPS) {
            p.x += velocity/yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
            p.y += velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
            p.theta += yaw_rate*delta_t;
        } else {
            p.x += velocity*cos(p.theta)*delta_t;
            p.y += velocity*sin(p.theta)*delta_t;
        }
        // add random Gaussian noise.
        p.x += x_dist(gen);
        p.y += y_dist(gen);
        p.theta += theta_dist(gen);
        
        particles[i] = p;
    }
    
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> preds, std::vector<LandmarkObs>& obs) {
    // Find the predicted measurement that is closest to each observed measurement
    for (int i=0; i<obs.size(); i++) {
        obs[i].id = 0;
        double min_dist = numeric_limits<double>::infinity();
        for (int p_id = 0; p_id<preds.size(); p_id++) {
            double dist_ = dist(obs[i].x, obs[i].y, preds[p_id].x, preds[p_id].y);
            // assign the observed measurement to this particular landmark.
            if (dist_ < min_dist) {
                obs[i].id = p_id;
                min_dist = dist_;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
