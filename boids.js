
// needed in constructor dx, boid_mass, boidspeed
function Flock(nboids,dx,boid_mass,boidspeed,vt_init,vr_init) {
  // An array for all the boids
  this.boid_set = []; // Initialize the array of boids
  this.dx = dx;  // pixels per unit length
  // this.maxforce = 20;
  // this.maxspeed = 2;
  this.boid_mass = boid_mass; // boid mass for flock
  this.boidspeed = boidspeed;
  this.vt_init = vt_init; // initial rotation of boids
  this.vr_init = vr_init; // initial random velocity dispersion of boids
  
  // no longer defined within flock
  // this.d_attract = d_attract; // distance for attraction forces between Boids
  // this.d_repel   = d_repel;  // distance for repulsion forces between Boids
  // this.d_align   = d_align;  // distance for alignment
  
  // this.attract_force = attract_force;  // attact force size
  // this.repel_force   = repel_force;    // repel force size
  // this.align_force_alpha   = align_force_alpha;    // steer, units of 1/time
  // this.propel_force_alpha  = propel_force_alpha;   // self propel, units 1/time

  // this.eta = eta; // for random motions
  // this.flip_v = flip_v;  // for flipping velocities
  
  // create nboids of Boids! constructor
  for (let i=0;i<nboids;i++){
    let b = new Boid(this.dx,this.boid_mass,this.vt_init, this.vr_init);
    this.boid_set.push(b);
  }
  
  this.display_flock = function(){  // display boids
    let n = this.boid_set.length;
    for (let i=0;i<n;i++){
      let b = this.boid_set[i];
      b.render();
    }
  }
  // call border function for each boid
  this.borders = function(){
    let n = this.boid_set.length;
    for (let i=0;i<n;i++){
      let b = this.boid_set[i];
      b.borders(); // rrap around border, periodic boundary
    }
  }
  
  // zero all the accelerations
  this.zeroaccel = function(){
    let n = this.boid_set.length;
    for(let k = 0;k<n;k++){
      this.boid_set[k].acceleration.mult(0);
    }
  }
  
  // Repel force depends on inverse distance -- all boid pairs
  // is a repel force if U_rep >0 otherewise is an attractive one
  // U_rep is now in units of energy
  // cuts off at distance d_rep
  // this version is 1/r force corresponding to potential
  // U = U_rep log r
  this.repel = function(d_rep,U_rep) {
   if (U_rep ==0){
      return;
   }
   let n = this.boid_set.length;
   for (let i = 0; i < n-1; i++) {
     let bi = this.boid_set[i];
     for (let j = i+1; j <n; j++) {
         let bj = this.boid_set[j];
         let dr = p5.Vector.sub(bi.position,bj.position);
         let r_len = dr.mag() + 0.01*d_rep;  // length of interboid distance
         // with softening
         if (r_len < d_rep){ // note cutoff distance here!!!!
           let drhat = dr.copy();
           drhat = drhat.normalize();  // unit vector for direction
           let Force = drhat.copy();
           let fac = -1.0*U_rep/r_len; // normalized here
           Force.mult(fac); 
           let ai = p5.Vector.mult(Force,-1/bi.m);
           let aj = p5.Vector.mult(Force, 1/bj.m);
           bi.acceleration.add(ai);
           bj.acceleration.add(aj);
        }
      // down side of this method is if we have many nearby particles, acceleration gets high
      }
    }
  }
  
  // Repel force depends on distance -- all boid pairs
  // is a repel force if U_rep >0 otherewise is an attractive one
  // U_rep is now in units of energy
  // cuts off at distance d_rep
  // this version is exponential from potential U = U_rep exp(-r/r_rep)
  this.repelexp = function(d_rep,U_rep) {
   if (U_rep ==0){
      return;
   }
   let n = this.boid_set.length;
   for (let i = 0; i < n-1; i++) {
     let bi = this.boid_set[i];
     for (let j = i+1; j <n; j++) {
         let bj = this.boid_set[j];
         let dr = p5.Vector.sub(bi.position,bj.position);
         let r_len = dr.mag() + 0.01*d_rep;  // length of interboid distance
         // with softening
         if (r_len < 2*d_rep){ // note cutoff distance here!!!!
           let drhat = dr.copy();
           drhat = drhat.normalize();  // unit vector for direction
           let Force = drhat.copy();
           let expfac = -1*U_rep*exp(-1.0*r_len/d_rep)/d_rep;
           Force.mult(expfac); 
           let ai = p5.Vector.mult(Force,-1/bi.m);
           let aj = p5.Vector.mult(Force, 1/bj.m);
           bi.acceleration.add(ai);
           bj.acceleration.add(aj);
        }
      }
    }
  }
  
  // exert a random turn on every particle each timestep or time called
  // eta in 1/time?
  this.random_eta = function(eta){
    if (eta == 0){
      return;
    }
    let n = this.boid_set.length;
    for (let i = 0; i < n; i++) {
        let bi = this.boid_set[i];
        let vi_perp = createVector(-1*bi.velocity.y,bi.velocity.x);
        vi_perp.normalize();
        vi_perp.mult(this.boidspeed);
        let noisesize = random(-eta, eta);
        let ai = p5.Vector.mult(vi_perp,noisesize/bi.m);  
        bi.acceleration.add(ai);     
    }
  }
  
  // flip velocity of a boid with probability flip_v in a single timestep
  this.flipv = function(flip_v){
    if (flip_v ==0){
      return;
    }
    let n = this.boid_set.length;
    for (let i = 0; i < n; i++) {
      let bi = this.boid_set[i];
      let ff = random(0,1.0);
      if (ff < flip_v){
        bi.velocity.mult(-1);
      }
    }
  }
  
    
  // Cohesion force, doing all pairs
  this.cohesion = function(attract_force,d_attract) {
     if (attract_force ==0){
       return;
     }
     let n = this.boid_set.length;
     for (let i = 0; i < n; i++) {
        let bi = this.boid_set[i];
        let pos_sum = createVector(0,0);
        let count = 0;
        for (let j = 0; j < n; j++) {
           let bj = this.boid_set[j];
           let dr = p5.Vector.sub(bi.position,bj.position);
           let r_len = dr.mag();  // length of interboid distance
           if ((r_len > 0) && (r_len < d_attract)){
              pos_sum.add(bj.position); // sum of positions
              count++;
           }
           if (count >0){
              pos_sum.div(count);   // is average of nearby boid positions!
              let Force = p5.Vector.sub(pos_sum, bi.position); // target direction
              Force.normalize();
              Force.mult(this.boidspeed); // boidspeed used here!
              Force.sub(bi.velocity);
              let ai = p5.Vector.mult(Force,attract_force/bi.m);  
              bi.acceleration.add(ai);
           }
        }
     }
  }
  
  // alignment try to stear toward mean velocity of nearby Boids
  // velocity dependent forces here!
  // nema gives nematic order version if it is 1
  this.align = function(nema,align_force_alpha,d_align) {  
     if (align_force_alpha==0){
       return;
     }
     let n = this.boid_set.length;
     // For every boid in the system, check if it's close to another
     for (let i = 0; i < n; i++) {
        let count = 0;
        let bi = this.boid_set[i];
        let vel_i = bi.velocity;
        // let steer = createVector(0,0);  // from average of nearest neighbor velocities
        let v_ave = createVector(0,0);
        for (let j = 0; j < n; j++) {
           let bj = this.boid_set[j];
           let dr = p5.Vector.dist(bi.position,bj.position);
           // If the distance is greater than 0 and less than an arbitrary amount (0 when you are yourself)
           if ((dr > 0) && (dr < d_align)) {
              let ifac = 1.0;
              let vel_j = bj.velocity.copy();
              if (nema ==1){ // nematic order
                let vidotvj = p5.Vector.dot(vel_i,vel_j);
                if (vidotvj < 0){
                    ifac = -1;
                }
              }
              vel_j.mult(ifac);
              v_ave.add(vel_j);  // sum of velocities of nearby
              count++; // Keep track of how many nearby
           }
        }
       
        if (count ==0){
           v_ave = bi.velocity.copy();
        }
       
        v_ave.normalize();
        v_ave.mult(this.boidspeed); // boidspeed used here!
        let steer = p5.Vector.sub(v_ave,bi.velocity);
        // desired is now the stear 
        // As long as the vector is greater than 0
        if (steer.mag() > 0) {
           steer.mult(align_force_alpha);
           // steer.limit(this.maxforce);
           bi.acceleration.add(steer);//  only on bi
        }
     }
  }
  
  // try to reach boid velocity, force is propto propel_force_alpha times 
  // difference in velocity from boidspeed
  // force is in direction of current velocity
  this.propel = function(propel_force_alpha){
    if  (propel_force_alpha ==0){
      return;
    }
    let n = this.boid_set.length;
    for (let i = 0; i < n; i++) {
       let bi = this.boid_set[i];
       let vi = bi.velocity.copy();
       let vmag = vi.mag();  // length
       let vhat = vi.copy(); // unit vector
       vhat.normalize();
       let ai = vhat.copy(); // direction same as velocity
       ai.mult((vmag-this.boidspeed)*propel_force_alpha*-1.0);
       bi.acceleration.add(ai);
    }
  }
  
  // integrate dt
  this.single_timestep = function(dt){
     // this.zeroaccel();
     // this.propel(); // propel boids
     // this.repel(); // repel repulsion
     // this.cohesion(); // attraction repulsion
     // this.align();
     // boid_node_interact(boid_set,node_set,force_amp,force_k,vforce_amp); 
     // apply interactions between boids and nodes
     let n = this.boid_set.length;
     for (let i = 0; i < n; i++) {
        bi = this.boid_set[i];
        let dv = p5.Vector.mult(bi.acceleration,dt); // h
        bi.velocity.add(dv); // 
        let dr = p5.Vector.mult(bi.velocity,dt);
        bi.position.add(dr); //
     }
  }
  
  // set all boid masses to this.boid_mass
  this.update_boidmass = function(){
     let n = this.boid_set.length;
     for (let i = 0; i < n; i++) {
       bi = this.boid_set[i];
       bi.m = this.boid_mass;
     }
  }
  
  // shift coordinates of all boids
  this.shift = function(centroid){
     let n = this.boid_set.length;
     for(let i = 0;i< n;i++){
        let boidi = this.boid_set[i];
        boidi.position.sub(centroid);
     }
     for(let i = 0;i< n;i++){ // move escaped boid back into rink
       let boidi = this.boid_set[i];
       if (abs(boidi.position.x) > 2){
          boidi.position.x = 0;
          boidi.velocity.x *=0.5;
       }
       if (abs(boidi.position.y) > 2){
          boidi.position.y = 0;
          boidi.velocity.y *=0.5;
       }
     }
     
  }
    
  // compute the mean v_theta  tangential velocity for all boids
  this.meanvtheta = function(){
    let n = this.boid_set.length;
    let sumv = 0.0;
    for(let i = 0;i< n;i++){
          let boidi = this.boid_set[i];
          let pos = boidi.position;
          let vel = boidi.velocity;
          let x = pos.x;
          let y = pos.y;
          let r = sqrt(x*x + y*y);
          let vtheta = (vel.y*x - vel.x*y)/r; // times theta hat
          sumv += vtheta;
     }
     let muv = sumv/n;
     return muv;
  }
  
  // compute the standard deviation from the velocity dispersion all boids
  this.stdv = function(){
    let n = this.boid_set.length;
    let sumvx = 0.0;
    let sumvy = 0.0;
    let sumv2 = 0.0
    for(let i = 0;i< n;i++){
          let boidi = this.boid_set[i];
          let vel = boidi.velocity;
          sumvx += vel.x;
          sumvy += vel.y;
          sumv2 += vel.x*vel.x + vel.y*vel.y;
     }
     let muvx = sumvx/n;
     let muvy = sumvy/n;
     let avev2 = sumv2/n;
     let sig2 = avev2 - muvx*muvx - muvy*muvy;
     return sqrt(sig2); // return standard deviation
  }
  
  
  // find the boid with the maximum number of neighbors 
  // within a distance rdist, return the number of nearest neighbors
  // for this boid
  this.maxnn = function(rdist){
    let n = this.boid_set.length;
    let imax = 0;
    let n_neighbors_max = 0;
    // let i_nn = 0;
    for(let i = 0;i< n;i++){
       let boidi = this.boid_set[i];
       let n_neighbors = 0;
       for(let j = 0;j< n;j++){
         let boidj = this.boid_set[j];
         let dr = p5.Vector.dist(boidi.position,boidj.position);
         if ((dr >0) && (dr < rdist)){
           n_neighbors++;
         }
       }
       if (n_neighbors > n_neighbors_max){
         n_neighbors_max = n_neighbors;
         // i_nn = i;
       }
    }
    return n_neighbors_max;
  }
  
  // reset all boids to initial conditions
  this.reset_boids = function(){
    let n = this.boid_set.length;
    for(let i = 0;i< n;i++){
      let boidi = this.boid_set[i];
      boidi.reset(this.vt_init,this.vr_init);
    }
  }
}

// Pos, Vel are the previously allocated vectors
// make Pos be a a randomly generated position
// evenly distributed in the  annulus enclosed by [rmin,rmax]
// also generate a velocity of rotation via vt, and with an
// an additional random component with size vt
// rewrites Pos,Vel from any old value if called as a reset 
function genposvel(rmin,rmax,vt,vr,Pos,Vel){
  let anglep = random(0,2*PI);
  let anglev = random(0,2*PI);
  let sr = random(rmin*rmin,rmax*rmax);
  let rw = sqrt(sr);
  Pos.x = rw*cos(anglep);
  Pos.y = rw*sin(anglep);
  Vel.x = -vt*sin(anglep);
  Vel.y =  vt*cos(anglep);
  Vel.x += vr*cos(anglev);
  Vel.y += vr*sin(anglev); 
}

// rmin,rmax generate in an annulus
const rmin_init = 0.5; // for initial conditions generation of boids, 
const rmax_init = 0.9; // also used for reseting initial conditions
// const vt_init = 0.0; // if you want rotation,  otherwise set to 0
// const vr_init = 0.1; // a random velocity component

function Boid(dxx,boid_mass,vt_init,vr_init) {
  // constructor and initialize with random position and velocity and rotation
  this.acceleration = createVector(0,0);
  this.velocity =  createVector(0,0);
  this.position =  createVector(0,0);
  genposvel(rmin_init,rmax_init,vt_init,vr_init,this.position,this.velocity);
  this.r = 4.0;  // for display and in pixels
  this.m = boid_mass;  
  this.dx = dxx; // scale multiply by this to get pixels
}

// reset initial condition for a single boid
// pos, vel choice is random
Boid.prototype.reset = function(vt_init,vr_init){
  genposvel(rmin_init,rmax_init,vt_init,vr_init,this.position,this.velocity);
}

// Boid display from p5.js examples
Boid.prototype.render = function() {
  // Draw a triangle rotated in the direction of velocity
  var theta = this.velocity.heading() + radians(90);
  push();
  fill(0,0,100);
  stroke(20);
  strokeWeight(1);
  translate(0,0);
  translate(width/2+this.position.x*this.dx,height/2+this.position.y*this.dx);
  rotate(theta);
  beginShape();
  vertex(0, -this.r*2);
  vertex(-this.r, this.r*2);
  vertex(this.r, this.r*2);
  endShape(CLOSE);
  pop();

}




// wrap boids around canvas border, but note that forces are not wrapped.
Boid.prototype.borders = function(){
  //print(this.dx);
  let wreal = (width + this.r)/this.dx ;
  let hreal = (height+ this.r)/this.dx ;
  if (this.position.x > 0.5*wreal){
    this.position.x -= wreal;
  }
  if (this.position.x < -0.5*hreal){
    this.position.x += wreal;
  }
  if (this.position.y > 0.5*wreal){
    this.position.y -= hreal;
  }
  if (this.position.y < -0.5*hreal){
    this.position.y += hreal;
  }
}
