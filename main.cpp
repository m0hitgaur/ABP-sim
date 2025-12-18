#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <format>
#include <thread>
#include <filesystem>
#include <iomanip>
#include "locker.h"
using namespace std;
namespace fs = filesystem;


bool create_directory(const string& path) {
    try {
        fs::create_directories(path);
        return fs::is_directory(path);
    } 
    catch (const fs::filesystem_error& e) {
        cerr << "Error creating directory '" << path << "': " << e.what() << '\n';
        return false;
    }
}
struct Particle {
    double x, y;
    double vx, vy;
    double ax, ay;
    double theta_new,theta;
    double x_new,y_new;
    double vx_new,vy_new;
    vector <int> neighbours;
};

class Simulation {
private:
    const int N;              // Number of particles
    const double Lx;          // Box size
    const double Ly;          //  Box size
    const double dt;          // Timestep
    const double v0;          // Magnitude of velocity
    const double Dr;          // Rotational diffusion coefficient
    const double Dt;          // Translational diffusion coefficient
    const double rc;          // Cutoff radius
    const double sigma;       // Particle radius
    const double k;           // Repulsion strength
    const int trial;          // Trial number    
    const int tmax;           // Total nujmber of time steps  
    string folder_path;  // Directory for saving data 
    vector<bool> time_record;  // Time points to record
    vector<int> times;         // Time points recorded 
    vector<double> order;     // Order parameter data
    vector<Particle> particles;
    mt19937 gen; 
    uniform_real_distribution<double> rng_uniform_symm{-1.0, 1.0};
    uniform_real_distribution<double> rng_uniform_one{0, 1};
    normal_distribution<double> rng_normal_unity{0.0, 1.0};

public:
    Simulation(int num_particles,double Dr_input,double Dt_input,double box_size_x,
               double box_size_y,double sigma_input,double k_input,double velo_mag,
               double timestep,int trial_input,int tmax_input,int seed)
               :N(num_particles), Dr(Dr_input),Lx(box_size_x),Ly(box_size_y), 
               dt(timestep), Dt(Dt_input),sigma(sigma_input), k(k_input),
               rc(3*sigma_input),v0(velo_mag),trial(trial_input),tmax(tmax_input){
        particles.resize(N);
        initialize_particles();
        initialize_time_to_record();
        folder_path="data/";
        create_directory(folder_path + "order_data");
        create_directory(folder_path+"config_data/trial_"+ to_string(trial)+"/");
        
    }

    vector<double> get_order_data(){
        return order;}
    vector<int> get_time_data(){
        return times;}

    void update_neigbours(){      
        for(Particle &p:particles){
            p.neighbours.clear();            
            for(int i=0;i<particles.size();i++){
                double r=rij(p,particles[i]);
                if(r<=2*rc && r>1e-5) {p.neighbours.push_back(i);}
            }
        }
    }

    void initialize_particles_old() {
        gen.seed(12345 + 10 * trial);
        double rho = N/(Lx*Ly);
        if(rho>=1){cout<<"Particles Intialize Failure ";}
        int i=0;
        for(int j=0;j<static_cast<int>(Lx) &&i<N;j++){
            for(int k=0;k<static_cast<int>(Ly)&&i<N;k++){
            
                particles[i].x = j ;
                particles[i].y = k;
                particles[i].x_new = j ;
                particles[i].y_new = k;                    
                
                double theta=M_PI*rng_uniform_symm(gen);
    
                particles[i].vx = v0*cos(theta);
                particles[i].vy = v0*sin(theta);
                particles[i].vx_new = v0*cos(theta);
                particles[i].vy_new = v0*sin(theta);
                particles[i].ax = 0;
                particles[i].ay = 0;
                i++;                
            } 
        }    
        update_neigbours();
        for(int t=0;t<100;t++){velocity_update();position_update();EndTimeStep();}
        update_neigbours();
    }                  
    void initialize_particles() {
        gen.seed(12345 + 10 * trial);
        double rho = N*M_PI*sigma*sigma/(Lx*Ly);
        if(rho>=1){cout<<"Density greater than 1  Packing Fraction = "<<rho<<endl<<flush;abort();}
        int grid_size = static_cast<int>(ceil(sqrt(N)));
        double spacing_x = Lx / grid_size;
        double spacing_y = Ly / grid_size;
        int i=0;
        for(int j=0;j<static_cast<int>(grid_size) &&i<N;j++){
            for(int k=0;k<static_cast<int>(grid_size)&&i<N;k++){

                particles[i].x = (j + sigma) * spacing_x;
                particles[i].y = (k + sigma) * spacing_y;
                particles[i].x_new = particles[i].x;
                particles[i].y_new = particles[i].y;
                
                double theta=M_PI*rng_uniform_symm(gen);
                particles[i].vx = 0;
                particles[i].vy = 0;
                particles[i].vx_new = 0;
                particles[i].vy_new = 0;
                
                particles[i].ax = 0;
                particles[i].ay = 0;
                
                particles[i].theta = theta;
                particles[i].theta_new = theta;
                i++;                
            } 
        }
        update_neigbours();
        //for(int t=0;t<1500;t++){velocity_update();position_update();EndTimeStep();}
        update_neigbours();
    }          
    double dot_product(double theta,double dx,double dy,double rij){
        return ( (cos(theta) * (dx)) + (sin(theta) * (dy)) )/(rij);
    }
    void pbc_position(Particle & p){
        // Periodic boundary conditions
        if (p.x_new < 0) p.x_new = fmod(p.x_new,Lx) +Lx;
        if (p.x_new > Lx) p.x_new = fmod(p.x_new,Lx);
        if (p.y_new < 0) p.y_new = fmod(p.y_new,Ly) +Ly;
        if (p.y_new > Ly) p.y_new = fmod(p.y_new,Ly);
    }
    double minimum_image(double dx,double L){
        if(dx>L/2) dx=dx-L;
        if(dx<-L/2) dx=dx+L ;  
        return dx;
    }
    void pbc_velo(Particle & p){
        double theta=atan2(p.vy_new,p.vx_new);
        if(theta<-M_PI)theta= fmod(theta , M_PI)+M_PI;
        else if(theta>M_PI)theta= fmod(theta , M_PI)-M_PI;
        p.vx_new=v0 *cos(theta);
        p.vy_new=v0 *sin(theta);
    }
    double lennard_jones_force(double r) {
        const double rmin = 1e-5;          
        if (r > rc || r < rmin) return 0.0;
        double sr  = sigma / r;
        double sr6 = pow(sr,6);
        double sr12 = sr6 * sr6;
        return 24.0 * (sr6 - 2.0 * sr12 ) / r;
    }
    double inter_particle_repulsive_force(double r) {
        if (r <= 2*sigma ) return k*(r-2*sigma);
        else return 0;
    }   
    double rij( Particle& p_i,  Particle& p_j) {
        double dx = p_i.x - p_j.x;
        double dy = p_i.y - p_j.y;
        dx=minimum_image(dx,Lx);
        dy=minimum_image(dy,Ly);
        double r=hypot(dx,dy);
        if (r < 1e-5) r = 1e-5;
        return r;
    }
    
    void compute_forces() {
        for (Particle &p : particles) {
            p.ax = 0;
            p.ay  = 0;
            }
        
        for (int i = 0; i < particles.size(); i++) {
            for (int j:particles[i].neighbours) {
                if(j<=i)continue; // to avoid double counting
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                dx=minimum_image(dx,Lx);
                dy=minimum_image(dy,Ly);
                double r = hypot(dx,dy);     
            
                double f = inter_particle_repulsive_force(r);
                double fx = f * (-dx / r);
                double fy = f * (-dy / r);
                
                particles[i].ax += fx;
                particles[i].ay += fy;
                
                particles[j].ax -= fx;
                particles[j].ay -= fy;
                
                
            }
        }
    }   
    void compute_forces_wo_neighlist() {
        for (Particle &p : particles) {
            p.ax = 0;
            p.ay  = 0;
            }
        for (int i = 0; i < particles.size(); i++) {
            for (int j=i+1;j<particles.size();j++) {
                //if(j<=i)continue; // to avoid double counting
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                dx=minimum_image(dx,Lx);
                dy=minimum_image(dy,Ly);
                double r = hypot(dx,dy);     
            
                double f = inter_particle_repulsive_force(r);
                double fx = f * (-dx / r);
                double fy = f * (-dy / r);
                
                particles[i].ax += fx;
                particles[i].ay += fy;
                
                particles[j].ax -= fx;
                particles[j].ay -= fy;
                
                
            }
        }
    }   
    void velocity_update(){
        compute_forces();
        for (Particle & p : particles) {    
            p.vx_new = p.ax + v0*cos(p.theta) + sqrt(2*Dt/dt)*rng_normal_unity(gen);
            p.vy_new = p.ay + v0*sin(p.theta) + sqrt(2*Dt/dt)*rng_normal_unity(gen); 
           
        }
    }
    void position_update(){
        for (Particle & p : particles) {    
            p.x_new += p.vx * dt ;
            p.y_new += p.vy * dt ;      
            p.theta_new+= sqrt(2*Dr*dt) * rng_normal_unity(gen);
            pbc_position(p);
        }
    }

    void EndTimeStep(){
       for (Particle &p : particles) {
            p.vx = p.vx_new;
            p.vy = p.vy_new;
            p.x = p.x_new;
            p.y = p.y_new;   
            p.theta = p.theta_new;         
            }         
    } 
    void integrate() {
        velocity_update();
        position_update(); 
        EndTimeStep();
    } 
    double velocity_order_parameter() {
        double sum_vx = 0, sum_vy = 0;
        
        for ( Particle & p : particles) {
            sum_vx += p.vx;         
            sum_vy += p.vy;    
        }
        
        double va = hypot(sum_vx,sum_vy) /(v0* N);
        return va;
    }   
    void save_snapshot(int step,int trial) {
        ofstream file(folder_path+"config_data/trial_"+to_string(trial)+"/config_" +to_string(step) + ".csv");
        
        file<<"x,y,vx,vy\n";
        for ( Particle & p : particles) {
            file << p.x << "," << p.y << ","
                 << p.vx << "," << p.vy << "\n";
               
        }
        
        file.close();
    }  


    ////////// order saving not working //////////////////////////
    void save_order_data(){
        ofstream order_file(folder_path+"order_data/order_parameter_"+to_string(trial)+"_.csv");
        order_file<<"va,t\n";
        for (int i=0;i<order.size();i++)order_file<<order[i]<<","<<times[i]<<"\n";      
        order_file.close();

    }
    ////////////////////////////////////////////////////////////
    void initialize_time_to_record(){
        for (int t = 0; t < tmax; t++) {
            bool should_record = false;
                if (t <= 10) should_record = true;
                if(t>=10 && t<100&& t%10==0) should_record = true;                   
                if (t<1000 && t >= 100 && t % 50 == 0) should_record = true;   
                if (t>=1000 && t % 100 == 0) should_record = true;                        
                time_record.push_back(should_record);
                }
        }

    void start_run() {
        ofstream f(folder_path+"parameters.csv");
        string head="N,Lx,Ly,Dr,v0,dt,Dt,maxiter,trial\n";
        f<< head;
        f<<N<<","<<Lx<<","<<Ly<<","<<Dr<<","<<v0
        <<","<<dt<<","<<Dt<<","<<tmax<<","<<trial; 
        f.close();
        double timestarted=static_cast <double>(time(NULL));
        int time_counter=0;
        int time_update_neighbours=static_cast<int>(1/(dt*4));
        for (int t = 0; t < tmax; t++) { 
            integrate();
            if(t%time_update_neighbours==0)update_neigbours();
            if (time_record[t]){
                order.push_back(velocity_order_parameter());    
                save_snapshot(t,trial);
                times.push_back(t);
                } 
            if (t % 500 == 0){print_progress(static_cast<double>(t)/static_cast<double>(tmax),timestarted);}
        } 
        print_progress(1.0,static_cast<double>(timestarted));
        save_order_data();
        cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << endl;
    }

};




int main() { 
    int N = 500;              // Number of particles
    double Lx = 24;         // Box size
    double Ly = 24;         // Box size
    double v0=1.0e0;          // Magnitude of velocity
    double dt = 1.0e-3;       // Timestep
    double Dr=1.0e-1 ;        // Rotational diffusion coefficient
    double Dt=1.0e0;         // Translational diffusion coefficient
    int tmax = 5.0e4;         // Maximum time
    int numberoftrials=1;     // Number of trials
    int trialstart=0;         // Starting trial number 
    double sigma=0.5;         // particle radius
    double k=120.0;            // repulsion strength
    int seed=12345;           // random seed
    time_t trial_time,start_time=time(NULL) , finish_time;
    for(int trial=trialstart;trial<numberoftrials;trial++){
        trial_time=time(NULL);
        cout<<"\n"<<"Trial number "<<trial<< " Out of "<<numberoftrials<<fixed<<setprecision(2)<<" || Packing Fraction : "<<((M_PI*N*(sigma)*(sigma))/(Lx*Ly))<<" | N = "<<N<<" || "<<endl ;
        Simulation sim(N,Dr,Dt,Lx,Ly,sigma,k,v0,dt,trial,tmax,seed);
        sim.start_run();
        cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds ";  
        
    }   
    
    cout<<"\n"<<"Total time elapsed : "<< time(NULL) - start_time <<" seconds ";
    
    return 0;
}