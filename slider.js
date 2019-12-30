

// create a slider and store a label and position for it so we
// can update text labeling during draw
// each slider should have a different index so they are 
// in a grid
function Tslider(index,label,min,max,value){
  const slider_xlen = 80;
  const slider_xlen_s = '80px';
  const x0 = 10;  
  const xsep = 10; // separation between sliders in x direction
  const y0 = 10;
  const ysep = 20;
  const nx = 6;   // number of sliders in x direction, for grid
  this.label = label;
  var xindex = index%nx;
  this.slider_x = xindex*(slider_xlen + xsep) + x0;
  this.slider_y = y0 + (int)(index/nx)*ysep;
  // this.slider_w = slider_len;
  // let max = value*2;  // note that maximum is set to twice value
  let step = (max-min)/50; // number of possible values
  
  // constructor
  this.slider = createSlider(min,max,value,step);
  this.slider.position(this.slider_x, this.slider_y);
  this.slider.style('width',  slider_xlen_s);
  // this.slider.style('background-color', 255); mojoy
  
}

// label the slider
Tslider.prototype.text = function(){
   fill(0); // black
   text(this.label,this.slider_x + 10,this.slider_y+20); // offset the label
}
