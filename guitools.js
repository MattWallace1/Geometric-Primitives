/**
 * Convert a comma-separated string of 3D vector
 * coordinates into a glMatrix.vec3 object
 * @param {string} s x,y,z string of vector coordinates
 * @return {glMatrix.vec3} Vector version of values
 */
function splitVecStr(s) {
  ret = [];
  s.split(",").forEach(function(x) {
      ret.push(parseFloat(x));
  });
  if (ret.length != 3) {
    alert("Must have 3 comma-separated coordinates in a vector!");
  }
  return glMatrix.vec3.fromValues(ret[0], ret[1], ret[2]);
}

/**
 * Convert an array into a comma-separated
 * list of values
 * @param {list} v List
 * @param {int} k Number of decimal places (default 2)
 */
function vecToStr(v, k) {
  if (k === undefined) {
      k = 2;
  }
  s = "";
  for (let i = 0; i < v.length; i++) {
      s += v[i].toFixed(k);
      if (i < v.length-1) {
          s += ",";
      }
  }
  return s;
}

PLOT_COLORS = {
  "C0":"#0066ff", "C1":"#ff9933", "C2":"#33cc33", "C3":"#cc00ff",
  "C4":"#ff3300", "C5":"#996633"
}; 


/**
 * Return Plotly plots of equal axes that contain a set of vectors
 * @param {list of glMatrix.vec3} vs
 * @return x, y, and z axes plots 
 */
function getAxesEqual(vs) {
  //Determine the axis ranges
  minval = 0;
  maxval = 0;
  for (let i = 0; i < vs.length; i++) {
      for (let j = 0; j < 3; j++) {
          if (vs[i][j] < minval){ minval = vs[i][j]; }
          if (vs[i][j] > maxval){ maxval = vs[i][j]; }
      }
  }
  return {
  x:{ x: [minval, maxval], y: [0, 0], z: [0, 0],
    mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'xaxis'
  },
  y:{ x: [0, 0], y: [minval, maxval], z: [0, 0],
    mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'yaxis'
  },
  z:{ x: [0, 0], y: [0, 0], z: [minval, maxval],
    mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'zaxis'
  }};
}

function getMousePos(canvas, evt) {
  let rect = canvas.getBoundingClientRect();
  return {
      X: evt.clientX - rect.left,
      Y: evt.clientY - rect.top
  };
}



/**
 * A class for managing the GUI for angle computation
 */
class AngleGUI {
  constructor() {
    this.setupMenu();
  }

  /**
   * Run the getAngle() function using the chosen vectors
   */
  getAngle() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    try{
      let res = getAngle(a, b, c);
      this.result = res.toFixed(3) + " degrees";
      if (this.plotVectors) {
        this.drawVectors();
      }
    }
    catch(err) {
      alert("There was an error with your getAngle code!");
      this.result = "error";
      console.log(err);
    }
  }

  /**
   * Setup a menu with inputs for 3 vectors, as well as
   * a button for computing and a text area for display
   */
  setupMenu() {
    let menu = new dat.GUI();
    this.menu = menu;
    this.a = "0,0,0";
    this.b = "0,0,0";
    this.c = "0,0,0";
    this.plotVectors = true;
    this.result = "";
    menu.add(this, "a");
    menu.add(this, "b");
    menu.add(this, "c");
    menu.add(this, "plotVectors");
    menu.add(this, "getAngle");
    menu.add(this, "result").listen();
  }

  drawVectors() {
    let plots = [];
    let labels = ["a", "b", "c"];
    let vecs = [splitVecStr(this.a), splitVecStr(this.b), splitVecStr(this.c)];
    // Plot endpoints
    for (let i = 0; i < 3; i++) {
      let v = vecs[i];
      let viz = { x: [v[0]], y: [v[1]], z: [v[2]],
        mode: 'markers+lines', line: {color: '#ffffff', width: 10},
        type: 'scatter3d', name: labels[i],
        marker: {color: PLOT_COLORS["C"+i], size: 10, symbol: 'circle'}
      };
      plots.push(viz);
    }

    // Plot vectors between endpoints
    let v1 = vecs[0];
    labels = ["aa", "ab", "ac"]
    for (let i = 1; i < 3; i++) {
      let v2 = vecs[i];
      let viz = { x: [v1[0], v2[0]], y: [v1[1], v2[1]], z: [v1[2], v2[2]],
        mode: 'lines', line: {color: '#000000', width: 10},
        type: 'scatter3d', name: labels[i],
      };
      plots.push(viz);
    }

    let axes = getAxesEqual([this.a, this.b, this.c]);
    plots.push(axes.x);
    plots.push(axes.y);
    plots.push(axes.z);
    let layout = {
      autosize: false, width: 500, height: 500,
      margin: { l: 0, r: 0, b: 0, t: 65 }
    };
    Plotly.newPlot('angleViz', plots, layout);  
  }
}


/**
 * A class for managing the GUI for parallel and perpendicular
 * projection
 */
class ProjGUI {
  constructor() {
    this.setupMenu();
  }

  /**
   * Run the getAngle() function using the chosen vectors
   */
  getProjections() {
    let u = splitVecStr(this.u);
    let v = splitVecStr(this.v);
    let parProj = vec3.create();
    let perpProj = vec3.create();
    // Call student parallel projection code
    try{
      parProj = projVector(u, v);
      this.parProj = vecToStr(parProj);
    }
    catch(err) {
      alert("There was an error with parallel projection!");
      this.parProj = "error";
      console.log(err);
    }

    // Call student perpendicular projection code
    try{
      perpProj = projPerpVector(u, v);
      this.perpProj = vecToStr(perpProj);
    }
    catch(err) {
      alert("There was an error with perpendicular projection!");
      this.perpProj = "error";
      console.log(err);
    }

    if (this.plotVectors) {
      this.drawProjections(parProj, perpProj);
    }
  }

  /**
   * Setup a menu with inputs for 3 vectors, as well as
   * a button for computing and a text area for display
   */
  setupMenu() {
    let menu = new dat.GUI();
    this.menu = menu;
    this.u = "0,0,0";
    this.v = "0,0,0";
    this.plotVectors = true;
    this.parProj = "";
    this.perpProj = "";
    menu.add(this, "u");
    menu.add(this, "v");
    menu.add(this, "plotVectors");
    menu.add(this, "getProjections");
    menu.add(this, "parProj").listen();
    menu.add(this, "perpProj").listen();
  }

  /**
   * Plot u, v, and the parallel/perpendicular projections using plot.ly
   * @param {glMatrix.vec3} parProj Parallel projection
   * @param {glMatrix.vec3} perpProj Perpendicular projection
   */
  drawProjections(parProj, perpProj) {
    let u = splitVecStr(this.u);
    let v = splitVecStr(this.v);
    var uviz = { x: [0, u[0]], y: [0, u[1]], z: [0, u[2]],
      mode: 'markers+lines', line: {color: PLOT_COLORS[0], width: 10},
      type: 'scatter3d', name: 'u',
      marker: {color: PLOT_COLORS[0], size: 4, symbol: 'circle'}
    };
    var vviz = { x: [0, v[0]], y: [0, v[1]], z: [0, v[2]],
      mode: 'markers+lines', line: {color: PLOT_COLORS[1], width: 10},
      type: 'scatter3d', name:'v',
      marker: {color: PLOT_COLORS[1], size: 4, symbol: 'circle'}
    };
    var projviz = { x: [0, parProj[0]], y: [0, parProj[1]], z: [0, parProj[2]],
      mode: 'markers+lines', line: {color: PLOT_COLORS[2], width: 10},
      type: 'scatter3d', name:'proj',
      marker: {color: PLOT_COLORS[2], size: 4, symbol: 'circle'}
    };
    var projperpviz = { x: [0, perpProj[0]], y: [0, perpProj[1]], z: [0, perpProj[2]],
      mode: 'markers+lines', line: {color: PLOT_COLORS[3], width: 10, arrowhead:7},
      type: 'scatter3d', name:'projperp',
      marker: {color: PLOT_COLORS[3], size: 4, symbol: 'circle'}
    };
    var axes = getAxesEqual([u, v, parProj, perpProj]);
    var data = [uviz, vviz, projviz, projperpviz, axes.x, axes.y, axes.z];
    var layout = {
      autosize: false, width: 500, height: 500,
      margin: { l: 0, r: 0, b: 0, t: 65 }
    };
    Plotly.newPlot('projViz', data, layout);
  }
}


/**
 * Return a set of plotly plot dictionary parameters
 * for plotting the vertices and edges of a triangle
 * @param {glMatrix.vec3} a First endpoint
 * @param {glMatrix.vec3} b Second endpoint
 * @param {glMatrix.vec3} c Third endpoint
 */
function getTrianglePlots(a, b, c) {
  let plots = [];
  let labels = ["a", "b", "c"];
  let vecs = [a, b, c];
  // Plot endpoints
  for (let i = 0; i < 3; i++) {
    let v = vecs[i];
    let viz = { x: [v[0]], y: [v[1]], z: [v[2]],
      mode: 'markers+lines', line: {color: '#ffffff', width: 10},
      type: 'scatter3d', name: labels[i],
      marker: {color: PLOT_COLORS["C"+i], size: 10, symbol: 'circle'}
    };
    plots.push(viz);
  }

  // Plot vectors between endpoints
  labels = ["ab", "bc", "ca"]
  for (let i = 0; i < 3; i++) {
    let v1 = vecs[i];
    let v2 = vecs[(i+1)%3];
    let viz = { x: [v1[0], v2[0]], y: [v1[1], v2[1]], z: [v1[2], v2[2]],
      mode: 'lines', line: {color: '#000000', width: 10},
      type: 'scatter3d', name: labels[i],
    };
    plots.push(viz);
  }
  return plots;
}



/**
 * A class for managing the GUI for area computation
 */
class AreaGUI {
  constructor() {
    this.setupMenu();
  }

  /**
   * Run the getTriangleArea() function using the chosen vectors
   */
  getTriangleArea() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    try{
      let res = getTriangleArea(a, b, c);
      this.result = res.toFixed(3);
      if (this.plotTriangle) {
        this.drawTriangle();
      }
    }
    catch(err) {
      alert("There was an error with your getTriangleArea code!");
      this.result = "error";
      console.log(err);
    }
  }

  /**
   * Setup a menu with inputs for 3 vectors, as well as
   * a button for computing and a text area for display
   */
  setupMenu() {
    let menu = new dat.GUI();
    this.menu = menu;
    this.a = "0,0,0";
    this.b = "0,0,0";
    this.c = "0,0,0";
    this.plotTriangle = true;
    this.result = "";
    menu.add(this, "a");
    menu.add(this, "b");
    menu.add(this, "c");
    menu.add(this, "plotTriangle");
    menu.add(this, "getTriangleArea");
    menu.add(this, "result").listen();
  }

  drawTriangle() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    let plots = getTrianglePlots(a, b, c);
    let axes = getAxesEqual([a, b, c]);
    plots.push(axes.x);
    plots.push(axes.y);
    plots.push(axes.z);
    let layout = {
      autosize: false, width: 500, height: 500,
      margin: { l: 0, r: 0, b: 0, t: 65 }
    };
    Plotly.newPlot('areaViz', plots, layout);  
  }
}




/**
 * A class for managing the GUI for area computation
 */
class AboveOrBelowGUI {
  constructor() {
    this.setupMenu();
  }

  /**
   * Run the getAboveOrBelow() function using the chosen vectors
   */
  getAboveOrBelow() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    let d = splitVecStr(this.d);
    try{
      this.result = getAboveOrBelow(a, b, c, d);
      if (this.drawVectors) {
        this.plotVectors();
      }
    }
    catch(err) {
      alert("There was an error with your getAboveOrBelow code!");
      this.result = "error";
      console.log(err);
    }
  }

  /**
   * Setup a menu with inputs for 4 vectors, as well as
   * a button for computing and a text area for display
   */
  setupMenu() {
    let menu = new dat.GUI();
    this.menu = menu;
    this.a = "0,0,0";
    this.b = "0,0,0";
    this.c = "0,0,0";
    this.d = "0,0,0";
    this.drawVectors = true;
    this.result = "";
    menu.add(this, "a");
    menu.add(this, "b");
    menu.add(this, "c");
    menu.add(this, "d");
    menu.add(this, "drawVectors");
    menu.add(this, "getAboveOrBelow");
    menu.add(this, "result").listen();
  }

  plotVectors() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    let d = splitVecStr(this.d);
    // Triangle plot
    let plots = getTrianglePlots(a, b, c);
    // d point plot
    let dviz = { x: [d[0]], y: [d[1]], z: [d[2]],
      mode: 'markers+lines', line: {color: '#ffffff', width: 10},
      type: 'scatter3d', name: "d",
      marker: {color: PLOT_COLORS["C3"], size: 10, symbol: 'circle'}
    };
    plots.push(dviz);
    let axes = getAxesEqual([a, b, c, d]);
    plots.push(axes.x);
    plots.push(axes.y);
    plots.push(axes.z);
    let layout = {
      autosize: false, width: 500, height: 500,
      margin: { l: 0, r: 0, b: 0, t: 65 }
    };
    Plotly.newPlot('aboveOrBelowViz', plots, layout);  
  }
}


class BarycentricGUI {
  constructor() {
    this.Ps = []; //Points [a, b, c] on the triangle
    this.p = null; // Point inside triangle
    this.coords = vec3.create();
    this.setupCanvas();
    this.selectTriangle();
  }

  /**
   * Setup the 2D canvas for selecting the triangle
   */
  setupCanvas() {
    let canvas = document.getElementById('barycanvas');
    let ctx = canvas.getContext("2d"); //For drawing
    ctx.font = "16px Arial";
    this.canvas = canvas;
    this.ctx = ctx;
    //Need this to disable that annoying menu that pops up on right click
    canvas.addEventListener("contextmenu", function(e){ e.stopPropagation(); e.preventDefault(); return false; }); 
    this.selectingTriangle = true;
    this.triangleButton = document.getElementById("selectingTriangle");
    this.pointButton = document.getElementById("selectingPoint");
    canvas.addEventListener("mousedown", this.selectVec.bind(this));
    canvas.addEventListener("touchstart", this.selectVec.bind(this)); //Works on mobile devices!
    this.repaint(); 
  }

  /**
   * Switch to selecting the point inside the triangle
   */
  selectPoint() {
      this.selectingTriangle = false;
      this.triangleButton.style.background = "#bfbfbf";
      this.pointButton.style.background = "#bb1111";
  }

  /**
   * Switch to selecting the points on the triangle
   */
  selectTriangle() {
      this.selectingTriangle = true;
      this.triangleButton.style.background = "#bb1111";
      this.pointButton.style.background = "#bfbfbf";
  }

  /**
   * Draw the triangle and the point inside it, as well
   * as the barycentric coordinates at each point
   */
  repaint() {
    let canvas = this.canvas;
    let ctx = this.ctx;
    let p = this.p;
    let Ps = this.Ps;
    let dW = 5;
    let W = canvas.width;
    let H = canvas.height;
    ctx.clearRect(0, 0, W, H);

    if (!(p === null)) {
        //Draw point
        ctx.fillStyle = PLOT_COLORS["C3"];
        ctx.fillRect(p[0]-dW, p[1]-dW, dW*2+1, dW*2+1);
    }
    
    //Draw triangle points
    for (let i = 0; i < Ps.length; i++) {
        ctx.fillStyle = PLOT_COLORS["C"+i];
        ctx.fillRect(Ps[i][0]-dW, Ps[i][1]-dW, dW*2+1, dW*2+1);
    }
    
    //Draw triangle edges
    ctx.fillStyle = "#000000";
    for (let i = 0; i < Ps.length; i++) {
        ctx.beginPath();
        ctx.moveTo(Ps[i][0], Ps[i][1]);
        ctx.lineTo(Ps[(i+1)%Ps.length][0], Ps[(i+1)%Ps.length][1]);
        ctx.stroke();        
    }

    // Draw coordinates next to the points
    if (Ps.length == 3) {
        this.coords.forEach(function(val, index) {
            ctx.fillText(val.toFixed(2), Ps[index][0] + 5, Ps[index][1]);
        });
    }
  }

  /**
   * Call the code and draw text where appropriate
   */
  finalizePointSelection() {
    let p = this.p;
    let Ps = this.Ps;
    if (Ps.length == 3 && !(p === null)) {
      try {
        this.coords = getBarycentricCoords(Ps[0], Ps[1], Ps[2], p);
      }
      catch(err) {
        alert("There was an error with your getBarycentricCoords code!");
        console.log(err);
      }
    }
    else {
        this.coords = this.coords.map(function(){return 0});
    }
    this.repaint();
  }

  selectVec(evt) {
      let mousePos = getMousePos(this.canvas, evt);
      let Ps = this.Ps;
      let p = this.p;
      let X = mousePos.X;
      let Y = mousePos.Y
      let clickType = "LEFT";
      evt.preventDefault();
      if (evt.which) {
          if (evt.which == 3) clickType = "RIGHT";
          if (evt.which == 2) clickType = "MIDDLE";
      }
      else if (evt.button) {
          if (evt.button == 2) clickType = "RIGHT";
          if (evt.button == 4) clickType = "MIDDLE";
      }
      
      if (this.selectingTriangle) {
          if (clickType == "LEFT") {
              //Add a point
              if (Ps.length < 3) {
                  Ps.push(vec3.fromValues(X, Y, 0));
              }
              else {
                  //If there's already a third point, simply replace it
                  Ps[2] = vec3.fromValues(X, Y, 0);
              }
          }
          else {
              //Remove point
              if (Ps.length > 0) {
                  Ps.pop();
              }
          }
          //Update text describing point coordinates
          for (let i = 0; i < 3; i++) {
              if (i < Ps.length) {
                  document.getElementById("Point"+i).innerHTML = "(" + Ps[i][0].toFixed(1) + "," + Ps[i][1].toFixed(1) + ")";    
              }
              else {
                  document.getElementById("Point"+i).innerHTML = "Not Selected";            
              }
          }
      }
      else {
          p = vec3.fromValues(X, Y, 0);
          this.p = p;
          document.getElementById("p").innerHTML = "(" + p[0].toFixed(1) + "," + p[1].toFixed(1) + ")";
      }
      this.finalizePointSelection();
  }
}

/**
 * A class for managing the GUI for ray triangle intersection
 */
class RayTriangleGUI {
  constructor() {
    this.setupMenu();
  }

  /**
   * Run the getAboveOrBelow() function using the chosen vectors
   */
  rayIntersectTriangle() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    let p0 = splitVecStr(this.p0);
    let v = splitVecStr(this.v);
    try{
      this.intersection = rayIntersectTriangle(p0, v, a, b, c);
      if (this.intersection.length > 0) {
        this.result = vecToStr(this.intersection[0]);
      }
      else {
        this.result = "[]";
      }
      if (this.drawVectors) {
        this.plotVectors();
      }
    }
    catch(err) {
      alert("There was an error with your getAboveOrBelow code!");
      this.result = "error";
      console.log(err);
    }
  }

  /**
   * Setup a menu with inputs for 4 vectors, as well as
   * a button for computing and a text area for display
   */
  setupMenu() {
    let menu = new dat.GUI();
    this.menu = menu;
    this.a = "0,0,0";
    this.b = "0,0,0";
    this.c = "0,0,0";
    this.p0 = "0,0,0";
    this.v = "0,0,0";
    this.drawVectors = true;
    this.result = "";
    let triMenu = menu.addFolder("Triangle");
    triMenu.add(this, "a");
    triMenu.add(this, "b");
    triMenu.add(this, "c");
    let rayMenu = menu.addFolder("Ray");
    rayMenu.add(this, "p0");
    rayMenu.add(this, "v");
    menu.add(this, "drawVectors");
    menu.add(this, "rayIntersectTriangle");
    menu.add(this, "result").listen();
  }

  plotVectors() {
    let a = splitVecStr(this.a);
    let b = splitVecStr(this.b);
    let c = splitVecStr(this.c);
    let p0 = splitVecStr(this.p0);
    let v = splitVecStr(this.v);
    // Triangle plot
    let plots = getTrianglePlots(a, b, c);
    // Vector plot
    var p0viz = { x: [p0[0]], y: [p0[1]], z: [p0[2]],
      mode: 'markers+lines', line: {color: '#ffffff', width: 10},
      type: 'scatter3d', name: 'p0',
      marker: {color: PLOT_COLORS["C3"], size: 10, symbol: 'circle'}
    };
    var vviz = { x: [p0[0], p0[0]+v[0]], y: [p0[1], p0[1]+v[1]], z: [p0[2], p0[2]+v[2]],
      mode: 'lines', line: {color: PLOT_COLORS["C4"], width: 10},
      type: 'scatter3d', name: 'v',
    };
    plots.push(p0viz);
    plots.push(vviz);
    // Intersection plot
    let res = this.intersection;
    if (res.length > 0) {
      var dviz = { x: [res[0][0]], y: [res[0][1]], z: [res[0][2]],
      mode: 'markers+lines', line: {color: "#ffffff", width: 10},
      type: 'scatter3d', name: 'd',
      marker: {color: PLOT_COLORS["C5"], size: 10, symbol: 'circle'}
      };

      var p0dviz = { x: [p0[0], res[0][0]], y: [p0[1], res[0][1]], z: [p0[2], res[0][2]],
      mode: 'lines', line: {color: PLOT_COLORS["C5"], width: 10},
      type: 'scatter3d', name: 'p0d',
      };
      plots.push(dviz);
      plots.push(p0dviz);
    }

    let axes = getAxesEqual([a, b, c, p0]);
    plots.push(axes.x);
    plots.push(axes.y);
    plots.push(axes.z);
    let layout = {
      autosize: false, width: 500, height: 500,
      margin: { l: 0, r: 0, b: 0, t: 65 }
    };
    Plotly.newPlot('rayTriViz', plots, layout);  
  }
}