#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <format>
#include <thread>
#include <filesystem>
#include <iomanip>
using namespace std;
namespace fs = filesystem;


bool create_directory(const string& path) {
    try {
        fs::create_directories(path);
        return fs::is_directory(path);
    } 
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directory '" << path << "': " << e.what() << '\n';
        return false;
    }
}
struct Particle {
    double x, y;
    double vx, vy;
    double ax, ay;
    double x_new,y_new;
    double vx_new,vy_new;
    vector <int> neighbours;
};

class Simulation_ABP {
private:
    const int N;              // Number of particles
    const double Lx;          // Box size
    const double Ly;          //  Box size
    const double dt;          // Timestep
    const double v0;          // Magnitude of velocity
    const double Temp;       // Temperture 
    const double D;          // Diffusion coefficient
    const double gamma;      // Friction coefficient
    const double rc;          // Cutoff radius
    const double sigma;       // Particle radius
    const double k;           // Repulsion strength
    const int trial;          // Trial number    
    vector<int> times;        // Time points recorded
    vector<double> order;     // Order parameter data
    vector<Particle> particles;
    mt19937 gen;
    uniform_real_distribution<double> rng_uniform_symm,rng_uniform_one;
    normal_distribution<double> rng_normal;    string folder_path;

public:
    Simulation_ABP(int num_particles, double Temp_strength,double D_input,double gamma_input,double box_size_x,
               double box_size_y,int sigma_input,double k_input,double velo_mag,double timestep,
               string folderpath,int trial_input,int seed)
             : N(num_particles), Lx(box_size_x),Ly(box_size_y), dt(timestep), Temp(Temp_strength),
               sigma(sigma_input), k(k_input),rc(3*sigma_input), rng_uniform_symm(-1.0, 1.0),rng_uniform_one(0, 1),
               gen(seed),rng_normal(0,1),v0(velo_mag),folder_path(folderpath),trial(trial_input),gamma(gamma_input),D(D_input){

        particles.resize(N);
        initialize_particles();
        
    }
    void save_order_data(){
        ofstream order_file(folder_path+"order_data/order_parameter_"+to_string(trial)+"_.csv");
        order_file<<"va,t\n";
        for (int i=0;i<order.size();i++)order_file<<order[i]<<","<<times[i]<<"\n";      
        order_file.close();

    }

    vector<double> get_order_data(){
        return order;}
    vector<int> get_time_data(){
        return times;}
    void update_neigbours(){
        
        for(Particle &p:particles){
            p.neighbours.clear();
            
            for(int i=0;i<N;i++){
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
        double rho = N/(Lx*Ly);
        if(rho>=1){cout<<"Particles Intialize Failure ";}
        int grid_size = static_cast<int>(ceil(sqrt(N)));
        double spacing_x = Lx / grid_size;
        double spacing_y = Ly / grid_size;
        int i=0;
        for(int j=0;j<static_cast<int>(grid_size) &&i<N;j++){
            for(int k=0;k<static_cast<int>(grid_size)&&i<N;k++){

                particles[i].x = (j + 0.5) * spacing_x;
                particles[i].y = (k + 0.5) * spacing_y;
                particles[i].x_new = particles[i].x;
                particles[i].y_new = particles[i].y;
                
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
        for(int t=0;t<1500;t++){velocity_update();position_update();EndTimeStep();}
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
    double minimum_image(double dx,double Lx){
        if(dx>Lx/2) dx=dx-Lx;
        if(dx<-Lx/2) dx=dx+Lx ;  
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
        if (r <= 2*sigma ) return k*(2*sigma-r);
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
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
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
            double theta=sqrt(2*D)*rng_normal(gen); 
            p.vx_new = p.ax/gamma + v0*cos(theta);
            p.vy_new = p.ay/gamma + v0*sin(theta); 
        }
    }
    void position_update(){
        for (Particle & p : particles) {    
            p.x_new += p.vx * dt ;
            p.y_new += p.vy * dt ;           
            pbc_position(p);
        }
    }

    void EndTimeStep(){
       for (Particle &p : particles) {
            p.vx = p.vx_new;
            p.vy = p.vy_new;
            p.x = p.x_new;
            p.y = p.y_new;            
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
        ofstream file_2(folder_path+"config_data/trial_"+to_string(trial)+"/force_" +to_string(step) + ".csv");
        file<<"x,y,vx,vy\n";
        for ( Particle & p : particles) {
            file << p.x << "," << p.y << ","
                 << p.vx << "," << p.vy << "\n";
            file_2<<p.ax<<","<<p.ay<<"\n";     
        }
        
        file.close();
    }   
    void start_run(int tmax,int trialstart,vector<bool>time_record) {
        ofstream f(folder_path+"parameters.csv");
        string head="N,Lx,Ly,v0,dt,eta,maxiter,trial\n";
        f<< head;
        f<<N<<","<<Lx<<","<<Ly<<","<<v0<<","<<dt<<","<<Temp<<","<<tmax<<","<<trial; 
        f.close();
        int time_counter=0;
        for (int t = 0; t < tmax; t++) { 
            compute_forces();
            integrate();
            if(t%250==0)update_neigbours();
            if (time_record[t]){
                order.push_back(velocity_order_parameter());    
                save_snapshot(t,trial);
                times.push_back(t);
                } 
            if (t % 500 == 0) cout  << t<<">>" << flush;
        } 
        cout<<tmax<<"\n"; 
        save_order_data();
        cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << endl;
    }

};

vector<bool> calculate_time_to_record(int tmax){
    vector<bool> time_record;
    for (int t = 0; t < tmax; t++) {
        bool should_record = false;
            if (t <= 10) should_record = true;
            if(t>=10 && t<100&& t%10==0) should_record = true;                   
            if (t<1000 && t >= 100 && t % 50 == 0) should_record = true;   
            if (t>=1000 && t % 100 == 0) should_record = true;                        
            time_record.push_back(should_record);
            }
    

    return time_record;
}
void save_order(vector<vector<double>>orderpara,vector<int>times,int trialstart,string folder_path){        
    ofstream order_file(folder_path+"order_parameter.csv");
    string a="";
    for (int i=0;i<orderpara.size();i++)a+="trial_"+to_string(i+trialstart)+",";      
    order_file<<a<<"t\n";
    for(int i=0;i<orderpara[0].size();i++){
        for(int j=0;j<orderpara.size();j++){
            order_file<<orderpara[j][i]<<",";
        }
        order_file<<times[i]<<"\n";
    }
    order_file.close();
    }

int main() {
    int N = 100;             // Number of particles
    double Lx = 20.0;          // Box size
    double Ly = 20.0;          // Box size
    double v0=50.0;            // Magnitude of velocity
    double D=0.1;            // Diffusion coefficient
    double gamma=1.0;        // Friction coefficient
    double dt = 1e-4;        // Timestep
    double Temp = 2.0;         // Temp 
    int tmax = 20000;        // Maximum time
    int numberoftrials=1;    // Number of trials
    int trialstart=0;        // Starting trial number 
    double sigma=1.0;          // particle radius
    double k=1;           // repulsion strength
    int seed=12345;          // random seed
    time_t trial_time,start_time=time(NULL) , finish_time;
    vector<bool> time_record=calculate_time_to_record(tmax); 
    vector <vector<double>> order_data(numberoftrials-trialstart);
    vector<int> times;
    string folder_path="data/";
    
    create_directory(folder_path + "order_data");
    for(int trial=trialstart;trial<numberoftrials;trial++){
        create_directory(folder_path+"config_data/trial_"+ to_string(trial)+"/");
        trial_time=time(NULL);
        cout<<"\n"<<"Trial number : "<<trial<< " Out of "<<numberoftrials<<fixed<<setprecision(2)<<" | D_r : "<<D<<" | Packing Fraction : "<<((N*M_PI*(sigma/2)*(sigma/2))/(Lx*Ly))<<" | N = "<<N<<" | "<<endl ;
        
        Simulation_ABP sim(N,Temp,D,gamma,Lx,Ly,sigma,k,v0,dt,folder_path,trial,seed);
        sim.start_run(tmax,trialstart,time_record);
        
        cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds ";  
        
        order_data[trial]=sim.get_order_data();       
        if(trial==numberoftrials-1)times=sim.get_time_data();
    }   
    
    save_order(order_data,times,trialstart,folder_path); 
    
    cout<<"\n"<<"Total time elapsed : "<< time(NULL) - start_time <<" seconds ";
    return 0;
}