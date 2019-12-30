// exponential short range forces between boids and nodes
// U propto force_amp*exp(-d/d_interact) where d = distance between
// force_amp is in units of force
function boid_node_interact(boid_set,node_set,
                             force_amp,d_interact,vforce_amp){
  let n_nod = node_set.length; // number of mass nodes in rink boundary
  let n_bd = boid_set.length;  // number of boids  
  
  for (let i=0;i<n_nod;i++){// loop over nodes in boundary
    let nodei = node_set[i];
    for(let j=0;j<n_bd;j++){ // loop over boids
      boidj = boid_set[j];
      let dr = p5.Vector.sub(boidj.position,nodei.position);  // vector between
      let d = dr.mag(); // distance between 
      //let v = dv.mag;
      if (d < 3*d_interact){
        let drhat = dr.copy();
        drhat.normalize(); //  is a unit vector
        let Force = drhat.copy(); // force in direction between
        Force.mult(-force_amp*exp(-d/d_interact));
        if (vforce_amp >0){
          let dv = p5.Vector.sub(boidj.velocity,nodei.velocity);
          let vForce = dv.copy();
          vForce.mult(vforce_amp);// damping force depends on velocity difference
          Force.add(vForce);
        }
        let a_boid = Force.copy();  // equal and opposite forces
        let a_node = Force.copy();
        a_boid.div(-boidj.m); // acceleration on boid
        a_node.div( nodei.m); // acceleration on node
        boidj.acceleration.add(a_boid); // apply to accelerations
        nodei.acceleration.add(a_node);
      }
    }
  }
  
}
