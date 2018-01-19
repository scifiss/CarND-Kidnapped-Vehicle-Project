/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *  Modified: Jan 1, 2018
 *      Author: Rebecca Gao
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

using namespace std;
default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;

	// Create a normal (Gaussian) distribution for x,y and theta from GPS data
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	double initWeight = 1.0;
	weights.resize(num_particles,initWeight);


	for (int i=0;i<num_particles;i++)
	{
	    Particle spart;
	    spart.id = i;
		 spart.x = dist_x(gen);
		 spart.y = dist_y(gen);
		 spart.theta = dist_theta(gen);
		 spart.weight = initWeight;
		 particles.push_back(spart);

	}

	is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Create a normal (Gaussian) distribution for x,y and theta for measurement
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	for (int i=0;i<num_particles;i++)
    {

        double theta = particles[i].theta;
        if (fabs(yaw_rate)>0.001)
        {
            particles[i].x += velocity/yaw_rate * (sin(theta + yaw_rate*delta_t)-sin(theta))+ dist_x(gen);
            particles[i].y += velocity/yaw_rate * (cos(theta) - cos(theta+ yaw_rate*delta_t)) + dist_y(gen);
            particles[i].theta += yaw_rate*delta_t+ dist_theta(gen);
        }
        else
        {
            particles[i].x += velocity*cos(theta)*delta_t+ dist_x(gen);
            particles[i].y += velocity*sin(theta)*delta_t+ dist_y(gen);
            particles[i].theta  = theta + dist_theta(gen);

        }

    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	double nearestdist, smalldist;
	int nearestid;
	for (int i=0;i<observations.size();i++)
    {
        nearestdist = 1e9;
        nearestid = -1;
        for (int j=0;j<predicted.size();j++)
        {
            smalldist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
            if (smalldist<nearestdist)
            {
                nearestdist = smalldist;
                nearestid = predicted[j].id;
            }

        }
        observations[i].id = nearestid;

    }


}

/**
compute 2D independent normal distribution pdf
*/
double normalIID2d_pdf(double x, double y, double mx, double my, double stdx, double stdy)
{
    double dx = (x-mx)/stdx;
    double dy = (y-my)/stdy;

    return exp( -( dx*dx + dy*dy)/2) /(2*M_PI*stdx*stdy);

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
	for (int i=0;i<num_particles;i++)
    {
        // observations from true car location for each particle
        vector<LandmarkObs> observ_map;
        //Particle spart=particles[i];

        // convert observations from each particle's coords to map's coords
        for (int j=0;j<observations.size();j++)
        {
            LandmarkObs obs;
            obs.id = -1;  // so far observations[j].id is not defined and not important
            obs.x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
            obs.y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
            observ_map.push_back(obs);
        }

        // get all landmarks that a particle can see
        vector<LandmarkObs> availmarks_map;
        for (int j=0;j<map_landmarks.landmark_list.size();j++)
        {

            //single_landmark_s landmark = map_landmarks.landmark_list[j];
            double dist_particle_landmark = dist(particles[i].x,particles[i].y,
                                                 map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f);
            if (dist_particle_landmark<=sensor_range)
            {
                LandmarkObs availableLandmark;
                availableLandmark.id = map_landmarks.landmark_list[j].id_i;
                availableLandmark.x = map_landmarks.landmark_list[j].x_f;
                availableLandmark.y = map_landmarks.landmark_list[j].y_f;
                availmarks_map.push_back(availableLandmark);
            }
        }
        // find in the observed marks's id in the map's landmark list
        dataAssociation(availmarks_map, observ_map);
//        if (observ_map.size()>availmarks_map.size())
//        {
//            cout << "sensed: " <<  observ_map.size() <<" avail: " << availmarks_map.size()<< endl;
//        }


        double pdftotal=1.0;
        if (availmarks_map.size()>0)
        {

            for (int j=0;j<availmarks_map.size();j++)
            {
                double pdfi;
                // in case different observations share one map landmark
                // detected observations are arithmetic averaged to get the middle position
                //  if 3 obs share one landmark
                // x = (x1+x2+x3)/3, y = (y1+y2+y3)/3
                int ndetected=0; // number of observed landmarks for each map landmark,
                double obs_x=0.0, obs_y=0.0;
                for (int k=0;k<observ_map.size();k++)
                {
                    if (availmarks_map[j].id == observ_map[k].id)
                    {
                        ndetected=ndetected+1;
                        //if (ndetected>0)
                        {
                        //    cout<<availmarks_map[j].id<<":"<<endl;
                        //    cout << availmarks_map[j].x<<","<<observ_map[k].x<<endl;
                        //    cout << availmarks_map[j].y<<","<<observ_map[k].y<<endl;
                        }
                        obs_x += observ_map[k].x;
                        obs_y += observ_map[k].y;

                    }
                }
                if (ndetected>0)
                {
                    if (ndetected>1)
                    {
                        //cout << "ndetected: " <<  ndetected << endl;
                        obs_x /= ndetected;
                        obs_y /= ndetected;
                        //cout << obs_x << ","<<obs_y<<endl;

                    }
                    pdfi = normalIID2d_pdf(obs_x, obs_y, availmarks_map[j].x, availmarks_map[j].y,
                                           std_landmark[0], std_landmark[1]);
//                    cout <<obs_x<<","<< obs_y<<","<<  availmarks_map[j].x<<","<<  availmarks_map[j].y<<
//                                           ":"<<pdfi<<endl;
                    pdftotal *= pdfi;
                }

            }
        }
        else // if the particle sees no landmark within sensor range, it's unlikely to be the car's position
            pdftotal = 1e-10;

        particles[i].weight = pdftotal;
        weights[i] = pdftotal;
        //cout<< "weight" << i <<":"<<pdftotal<<endl;
    }

}
/**
 * resample Resamples from the updated set of particles to form
 *   the new set of particles.
 */
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	discrete_distribution<int> weighti(weights.begin(),weights.end());
	vector<Particle> resampledParticles ;
	for (int i=0;i<num_particles;i++)
    {
//        Particle p
//        {
//            i,
//            particles[i].x,
//            particles[i].y,
//            particles[i].theta,
//
//
//        }
        resampledParticles.push_back(particles[weighti(gen)]);
    }

    particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
