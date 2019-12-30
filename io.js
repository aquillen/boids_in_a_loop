
// uses globlas nnodes, nboids, radfac, dt, big_radius, canvas_size
// U_repel, d_repel, d_align, align_force_alpha, force_amp, d_interact
function print_params_tofile(flock,mass_spring_system,filename){
  let writer = createWriter(filename);
// creates a file in downloads, comma separated format
  writer.print(['boidspeed', flock.boidspeed]);
  writer.print(['big_radius', big_radius]);
  let MM_boids = flock.boid_mass*nboids;
  let MM_nodes = mass_spring_system.node_mass*nnodes;
  writer.print(['M_boids', MM_boids]);
  writer.print(['canvas_size', canvas_size]);
  writer.print(['rad_fac', rad_fac]);
  writer.print(['dt', dt]);
  writer.print(['nnodes', nnodes]);
  writer.print(['nboids', nboids]);
  writer.print(['node_mass', mass_spring_system.node_mass]);
  writer.print(['boid_mass', flock.boid_mass]);
  // writer.print(['ds', ds]);
  writer.print(['nema', nema]);
  writer.print(['align_force_alpha', align_force_alpha]);
  writer.print(['d_align', d_align]);
  writer.print(['mass_spring_system.gamma_node',mass_spring_system.gamma_node]);
  writer.print(['mass_spring_system.ks', mass_spring_system.ks]);
  // writer.print(['U_attract', U_attract/boid_mass]);
  // writer.print(['d_attract', d_attract]);
  writer.print(['U_repel_v2', U_repel/flock.boid_mass]);
  writer.print(['d_repel', d_repel]);
  writer.print(['mass_spring_system.b_alpha', mass_spring_system.b_alpha]);
  // writer.print(['eta', eta]);
  // writer.print(['flip_v',flip_v]);
  writer.print(['force_amp',force_amp]);
  writer.print(['d_interact',d_interact]);
  writer.print(['Mass_ratio_boundary_boid',MM_nodes/MM_boids]);
  let epsilon_ks = nboids*flock.boid_mass*flock.boidspeed*flock.boidspeed;
  epsilon_ks /= (big_radius*mass_spring_system.ks*mass_spring_system.ds);
  writer.print(['epsilon_ks',epsilon_ks]);
  
  let pmax = flock.maxnn(0.5)/nboids;
  writer.print(['percent_max_nn',pmax]);
  let vstd = flock.stdv();
  writer.print(['vel_std',vstd]);
  let mvt_boid = flock.meanvtheta();
  writer.print(['mean_vt_boid',mvt_boid]);
  let mvt_node = mass_spring_system.meanvtheta();
  writer.print(['mean_vt_node',mvt_node]);
  
  // compute fourier coeffs
  let A2 = mass_spring_system.A_m(2);
  let A3 = mass_spring_system.A_m(3);
  let A4 = mass_spring_system.A_m(4);
  let A5 = mass_spring_system.A_m(5);
  let A6 = mass_spring_system.A_m(6);
  let A7 = mass_spring_system.A_m(7);
  let A47 = sqrt(A4*A4 + A5*A5 + A6*A6 + A7*A7);
  writer.print(['A2',A2]);
  writer.print(['A3',A2]);
  writer.print(['A47',A47]);
  
  let Pi_alpha = mass_spring_system.b_alpha/(MM_boids*flock.boidspeed*flock.boidspeed*big_radius);
  let Pi_lambda = sqrt(Pi_alpha*(MM_boids/MM_nodes))*pow(2.0*PI,1.5);
  let Pi_A = nboids*(d_repel/big_radius) * (d_repel/big_radius)/PI;
  let Pi_gas = align_force_alpha*d_align/flock.boidspeed;
  let Pi_push = (MM_boids/MM_nodes)*(align_force_alpha/mass_spring_system.gamma_node);
  let Pi_U = sqrt((MM_boids/MM_nodes)*(U_repel/(flock.boid_mass*flock.boidspeed*flock.boidspeed))*2*PI);
  
  writer.print(['Pi_alpha',Pi_alpha]);
  writer.print(['Pi_lambda',Pi_lambda]);
  writer.print(['Pi_A',Pi_A]);
  writer.print(['Pi_gas',Pi_gas]);
  writer.print(['Pi_push',Pi_push]);
  writer.print(['Pi_U',Pi_U]);
  
  writer.close();  // close file
  // saveCanvas('myCanvas', 'png'); // give it a try!
}

// print out boundary node positions and velocities
function print_bfile(mass_spring_system,filename){
  let writer = createWriter(filename);
  // creates a file in downloads, comma separated format

  let n = mass_spring_system.node_set.length;
  for (let i=0;i<n;i++){
      let nodei = mass_spring_system.node_set[i];
      let posi = nodei.position;
      let veli = nodei.velocity;
      writer.print([i,nodei.m,posi.x,posi.y,veli.x,veli.y]);
   }
   writer.close();  // close file
}
