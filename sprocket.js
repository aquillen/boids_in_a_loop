
var mass_spring_system; // boundary
var flock; // flock o boids
var sim_time = 0.0;  // simulation time 

var ofilename='Boids_sprocket.txt'; // for output table, will go into Downloads
var canvasname = 'Boids_sprocket'; // for snaps, png images

// units: ring radius and boid velocity
const boidspeed = 1;  // speed of self propulsion
const big_radius = 1.0;  // units of rink size 

// canvas and rink (boundary)
var canvas_size = 700; // pixels
const rad_fac = 0.60; // sets radius of rink boundary w.r.t to half of canvas_size
var rad_rink_px = rad_fac*canvas_size/2;  // in pixels
const dx = rad_rink_px/big_radius;  // grid spacing this is pixel length
// dx is only used for display

// timestep
const dt = 5e-3*big_radius/boidspeed; // timestep 

// numbers of nodes/boids and mass ratio
const nnodes  = 150; // numbers of mass nodes in rink
const M_boundary  = 25.0; // total boundary mass
const node_mass = M_boundary/nnodes; // mass of nodes
const nboids = 400;  // number of boids
const M_boids = 1.0; // total boid mass
const boid_mass = M_boids/nboids; // boid mass

// bending force 
const b_alpha = 1e-3;  // bending force units is energy times distance

// damping on boundary nodes //
const gamma_node = 0.1;  // damping on all boundary nodes 
// (force depends on velocity, units 1/time because gives acceleration

////////////////////////
//  boid forces 

var d_align = 0.2; // scale for alignment
var align_force_alpha = 3.0; // units of 1/time!

var d_repel = 0.1; // scale for repel force
var U_repel = 0.1*boid_mass;  // Urepel: units is energy

// boid/boundary interaction /////////////////
const force_amp = 1.5; // for interactions between boids and nodes
                      // is exponential -- this is units of force
const d_interact = 0.02;  // scale for interactions between boids and nodes
const vforce_amp = 0.00;  // damping between nodes and boids

//springs and boundary
var ks = 2e4*node_mass; // spring constant, units force per unit length
const gammas = 0.0;   // spring damping parameter

// initial conditions
const vt_init = 0.8;  // initial rotation  of boids
const vr_init = 0.1;  // random velocity for initial conditions

///////////////////// below not to be set or not used

const ten_speed = 0e0; // tension bending force, units speed^2, not used
var nema = 0; // if 0 self-propelled,  if 1 then nematic
var eta = 0e-1;  // noise
var flip_v = 0e-3;  // probability of flipping v direction in a timestep

////////////////used!
// display
const ndt=3;  // number of time steps per display update

var onscreen = 0; /// show some numbers on canvas
var usesliders = 0; // use sliders or not

var nbcounter = 0; // number of boundary files printed out
const pbinfo = 1;  // print boundary info or not

////////////////////////////////////////
// sliders 
var repel_force_Tslider; 
var d_repel_Tslider;
var align_force_Tslider;
var d_align_Tslider;
var gamma_node_Tslider;
var ks_Tslider;
var attract_force_Tslider; 
var d_attract_Tslider;
var b_alpha_Tslider;
var tension_Tslider;
var eta_Tslider;
var flip_Tslider;

function setup() {
  createCanvas(canvas_size, canvas_size); 
  // let xwidth = width*dx;     // grid size set here!
  // let yheight = height*dx;
  
  if (nema==0){
    eta = 0;
    flip_v = 0;
  }
  sim_time = 0.0; // simulation run time

  // set up masses/springs
  mass_spring_system = new Mass_spring_system(nnodes,big_radius,node_mass,dx);  
  mass_spring_system.display_springs();  // display springs
  mass_spring_system.display_nodes();   // display nodes
  
  // set up flock  of boids
  flock = new Flock(nboids,dx,boid_mass,boidspeed,vt_init,vr_init);  
  flock.display_flock(); // display flock of boids
  
  // set up some sliders
  if (usesliders==1){
    repel_force_Tslider = new Tslider(0,'repel',0,2*U_repel,U_repel); 
    d_repel_Tslider = new Tslider(1, 'd_repel',0,0.5,d_repel);
    align_force_Tslider = new Tslider(2, 'align',0,2*align_force_alpha,align_force_alpha);
    d_align_Tslider = new Tslider(3, 'd_align',0,0.5,d_align);
    gamma_node_Tslider = new Tslider(4, 'nodedamp',0,2*gamma_node,gamma_node);
    ks_Tslider = new Tslider(5, 'springk',ks/8,2*ks,ks);
    b_alpha_Tslider = new Tslider(6, 'alpha',0,5*b_alpha,b_alpha);
  
    //  attract_force_Tslider = new Tslider(6,'attract',0,2*U_attract,U_attract); 
    //  d_attract_Tslider = new Tslider(7, 'd_attract',0,2*d_attract,d_attract);
    //  tension_Tslider = new Tslider(9, 'tension',0,2*ten_speed,ten_speed);
    //  eta_Tslider = new Tslider(10, 'eta',0,2*eta,eta);
    //  flip_Tslider = new Tslider(11, 'flipv',0,2*flip_v,flip_v);
  }

 
}

var dcount=0; // used for centroiding display

function draw() {

   let boid_set = flock.boid_set;
   let node_set = mass_spring_system.node_set;
   background(240);
    
   mass_spring_system.display_springs(); // display springs
   mass_spring_system.display_nodes();  // display nodes
   flock.display_flock();  // display boid flock
    
   for (let k=0;k<ndt;k++){ // numbers of timesteps per display
     
     // zero accelerations
      mass_spring_system.zeroaccel();
      flock.zeroaccel();   
    
      // compute accelerations on mass/spring system
      mass_spring_system.compute_accel();
      
      // compute accelerations on boid flock
      flock.align(nema,align_force_alpha,d_align);
      flock.repelexp(d_repel,U_repel);
      // flock.propel(propel_force_alpha);
      // flock.cohesion(attract_force,d_attract);
      // flock.repelexp(d_attract,-1*U_attract);
      // flock.random_eta(eta);
      // flock.flipv(flip_v);
      
      // compute interactions boids/masses
      boid_node_interact(boid_set,node_set,force_amp,d_interact,vforce_amp);
    
      // update, do timestep!
      mass_spring_system.single_timestep(dt);
      flock.single_timestep(dt);
      sim_time += dt; // keep track of simulation time
  
      if (nnodes<2){
        flock.borders(); // if there aren't any mass nodes
      }
      else{ // centroid display
        dcount++;
        if ((dcount%1)==0){  // shift centroid
          let centroid = mass_spring_system.centroid();
          mass_spring_system.shift(centroid);
          flock.shift(centroid);
        }
      }
      
     // print boundary information into a file
     if (pbinfo==1){
       if (int(sim_time/dt) % 200 == 0){
         let bfile  = 'b' ;
         if  (nbcounter < 10){
           bfile =  bfile + '0'
         }
         if (nbcounter < 100){
           bfile =  bfile + '0'
         }
         bfile = bfile + str(nbcounter)  + '.txt'
         print_bfile(mass_spring_system,bfile)
         nbcounter +=1 ;
       }
     }
     
   }
   
   // use slider info! update flock and springsystem params
   if (usesliders==1){
     U_repel = repel_force_Tslider.slider.value();
     align_force_alpha = align_force_Tslider.slider.value();
     d_repel = d_repel_Tslider.slider.value();
     d_align = d_align_Tslider.slider.value();
     mass_spring_system.gamma_node = gamma_node_Tslider.slider.value();
     mass_spring_system.ks = ks_Tslider.slider.value();
     mass_spring_system.updateks(); // propagate to all springs
     mass_spring_system.b_alpha = b_alpha_Tslider.slider.value();
     // U_attract = attract_force_Tslider.slider.value();
     // d_attract = d_attract_Tslider.slider.value();
     // mass_spring_system.ten_speed = tension_Tslider.slider.value();
     // eta = eta_Tslider.slider.value();
     // flip_v = flip_Tslider.slider.value();
   
     repel_force_Tslider.text();
     align_force_Tslider.text();
     d_repel_Tslider.text();
     d_align_Tslider.text();
     gamma_node_Tslider.text();
     ks_Tslider.text();
     b_alpha_Tslider.text();
     // attract_force_Tslider.text();
     // d_attract_Tslider.text();
     // tension_Tslider.text();
     // eta_Tslider.text();
     // flip_Tslider.text();
   }
   if (onscreen==1){ // print stuff on screen
     textAlign(LEFT);
     let pmax = trunc(flock.maxnn(0.5)/nboids);
     text('pmax:' + str(pmax),10,70);
     let vstd = flock.stdv();
     text('vstd:' + str(trunc(vstd)),10,70+20);
     let mvt_boid = flock.meanvtheta();
     text('mvt_boid:' + str(trunc(mvt_boid)),10, 70+40);
     let mvt_node = mass_spring_system.meanvtheta();
     text('mvt_node:' + str(trunc(mvt_node)),10, 70+60);
     let A2 = mass_spring_system.A_m(2);
     let A3 = mass_spring_system.A_m(3);
     let A4 = mass_spring_system.A_m(4);
     let A5 = mass_spring_system.A_m(5);
     let A6 = mass_spring_system.A_m(6);
     let A7 = mass_spring_system.A_m(7);
     let A47 = sqrt(A4*A4 + A5*A5 + A6*A6 + A7*A7);
     text('A2:' +str(trunc(A2 )),10, 70+80);
     text('A47:'+str(trunc(A47)),10, 70+100);
   }
   text('t:' + str(trunc(sim_time)),10, canvas_size-20);
   
// how to reset simulation
//     flock.reset_boids();
//     mass_spring_system.reset();
//  sim_time = 0.0   


}


var stopstart=0;
function keyTyped() {
  if (key === 'w') {    
    print('hello');
    print_params_tofile(flock,mass_spring_system,ofilename);
    saveCanvas(canvasname, 'png'); 
  } 
  if (key === 'q') {  
    if (stopstart==0){
     noLoop(); // stop simulation
     stopstart=1;
    }
    else {
      loop();
      stopstart=0;
    }
  } 
  if (key === 's'){ // toggle  screen output
  
    onscreen += 1; /// show some numbers on canvas
    onscreen = onscreen%2;
  
    
  }
  
}

// just for printing numbers nicely 
function trunc(x){
  let z = 0.00000 + int(x*1000)/1000;
  return z;
}

function trunc_n(x,n){
  let p10 = pow(10.0,n);
  let z = 0.00000 + int(x*p10)/p10;
  return z;
}
