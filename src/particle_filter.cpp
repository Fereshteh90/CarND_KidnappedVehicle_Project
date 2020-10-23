/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  //default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  //creates a normal (Gaussian) distribution for x and y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta,std_theta);
  
   //particles.resize(num_particles); 

   
  // Sample from these normal distributions
  for (int i = 0; i < num_particles; ++i) {
    
    Particle ptc;
    ptc.id = i;
    ptc.x = dist_x(gen);
    ptc.y = dist_y(gen);
    ptc.theta = dist_theta(gen);
    ptc.weight = 1.0;
    
    particles.push_back(ptc);
    weights.push_back(1.0);  
  }
  
  is_initialized = true;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  //creates a normal (Gaussian) distribution
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0,std_pos[2]);
  
  //default_random_engine gen;
  for (int i = 0; i < num_particles; ++i)
  {
    // when the yaw rate is close to zero:
    if (fabs(yaw_rate) < 0.00001)
    {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    // when the yaw rate is not zero:
    else
    {
      particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
    // adding noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    
  }
    
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (unsigned int i = 0; i < observations.size(); i++)
  {
    // initial the closest distance
    double close_dist = numeric_limits<double>::max();
    int land_close_id = 0;
    
     for (unsigned int j = 0; j < predicted.size(); j++)
     {
       // calculate distance between observation and prediction 
       double estimate_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
       
       if (estimate_dist < close_dist)
       {
         close_dist = estimate_dist;
         observations[i].id = land_close_id;            
       } 
       land_close_id ++;
     }
 
  }//observations[i].id = land_close_id;

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  for(int i=0; i < num_particles; i++)
  {
    
    // step 1: transformed observation in the car coordinate to map coordinate
     vector<LandmarkObs> transformed_obs;
     transformed_obs.clear();
    
     for (int j=0; j<observations.size(); j++)
     {
      //Homogenous Transformation
     double obs_mapX = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) *                              observations[j].y;
     double obs_mapY = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) *                              observations[j].y;  
     transformed_obs.push_back(LandmarkObs{observations[j].id, obs_mapX, obs_mapY}); 
     }
    
    
  // step 2: landmarks in particles range 
    vector<LandmarkObs> landmarks_choice;
    for (unsigned j=0; j < map_landmarks.landmark_list.size(); j++)
    {
      double landmarks_x = map_landmarks.landmark_list[j].x_f;
      double landmarks_y = map_landmarks.landmark_list[j].y_f;
      int landmarks_id = map_landmarks.landmark_list[j].id_i; 
      
      if (dist(particles[i].x, particles[i].y, landmarks_x, landmarks_y) <= sensor_range)
      {
       landmarks_choice.push_back(LandmarkObs{landmarks_id, landmarks_x, landmarks_y});  
      }
   
    }
    
     if (landmarks_choice.size() == 0)
     {
      
      continue; 
     }
  
    
    // step 3: data association 
    dataAssociation(landmarks_choice, transformed_obs);
    
    // step 4: update weight 
    particles[i].weight = 1.0;
    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double cov_x = std_x * std_x;
    double cov_y = std_y * std_y;
    // calculate normalization term
    double gauss_norm = 1 / (2 * M_PI * std_x * std_y);
    
    for (unsigned int j=0; j<transformed_obs.size(); j++)
    {
      double x_obs = transformed_obs[j].x;
      double y_obs = transformed_obs[j].y;
      int id_land = transformed_obs[j].id;
      //double mu_x;
      //double mu_y;
      //for(size_t k=0; k<landmarks_choice.size(); k++)
      //{
        //if(landmarks_choice[k].id == transformed_obs[j].id)
        //{
         //mu_x = landmarks_choice[k].x;
         //mu_y = landmarks_choice[k].y;
        //}
      //} 
      double mu_x = landmarks_choice[id_land].x;
      double mu_y = landmarks_choice[id_land].y;
      double power_x = pow((x_obs - mu_x),2);
      double power_y = pow((y_obs - mu_y),2);
      double exponent = (power_x/(2 * cov_x)) + (power_y/(2 * cov_y));
      //double exponent = (pow(x_obs - mu_x, 2) / (2 * cov_x)) + (pow(y_obs - mu_y, 2) / (2 * cov_y));
      particles[i].weight *= gauss_norm * exp(-exponent);
      //particles[i].weight = particles[i].weight * weight_est ;
      cout << "weight: " << particles[i].weight << endl;
    }
    
    weights[i] = particles[i].weight;
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> resample_particles;
  uniform_real_distribution<> uni_dist(0, 1);
  int index = rand() % num_particles;
  double beta = 0.0;
  const double mw = *max_element(weights.cbegin(), weights.cend());
  
  for (int i=0; i < num_particles; i++)
  {
    beta += uni_dist(gen) * 2.0 * mw;
    while (beta > weights[index]) 
    {
      beta -= weights[index];
      index = (index + 1) % num_particles;
     }
      
    resample_particles.push_back(particles[index]);
  }
  
   particles = resample_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}