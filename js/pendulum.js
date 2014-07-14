var GRAVITY = 9.81; // meters / second^2
var DENSITY = 999.97; // kilograms / meter^3

var INTEGRATION_METHODS = {
  FORWARD_EULER: { names: ['Forward Euler', 'forward'] },
  SEMI_IMPLICIT_EULER: { names: ['Semi-implicit Euler', 'semi'] },
  RUNGE_KUTTA: { names: ['Runge-Kutta', 'runge-kutta'] },
};

// lengths are in meters
// masses are in kilograms
// thetas are in radians
// omegas are radians / second
var Pendulum = function(options) {
  var allowed_options = [
    'lengths', 'masses', 'thetas', 'omegas', 'integration_method'
  ];
  _.extend(this, _.pick(options, allowed_options));
  if (_.isUndefined(this.integration_method)) {
    this.integration_method = INTEGRATION_METHODS.SEMI_IMPLICIT_EULER;
  }

  if (this.lengths.length !== this.masses.length
      || this.masses.length !== this.thetas.length
      || this.thetas.length !== this.omegas.length) {
    throw "All parameters must be arrays of the same length.";
  }

  this.dimension = this.lengths.length;

  this.radii = [];
  _.each(this.masses, function(m) {
    this.radii.push(Math.pow( (3 * m) / (4 * Math.PI * DENSITY), 1 / 3));
  }, this);

  this.total_length = 0;
  _.each(this.lengths, function(len) { this.total_length += len; }, this);
  this.total_length += this.radii[this.radii.length - 1];

  this.time = 0;
};

Pendulum.prototype.getCartesianPositions = function() {
  var pos = [];
  var sx = 0;
  var sy = 0;
  for (var i = 0; i < this.dimension; i++) {
    sx += this.lengths[i] * Math.sin(this.thetas[i]);
    sy -= this.lengths[i] * Math.cos(this.thetas[i]);
    pos.push([sx, sy]);
  }
  return pos;
};

Pendulum.prototype.getCartesianVelocities = function() {
  var vel = [];
  var vx = 0;
  var vy = 0;
  for (var i = 0; i < this.dimension; i++) {
    vx += this.lengths[i] * Math.cos(this.thetas[i]) * this.omegas[i];
    vy += this.lengths[i] * Math.sin(this.thetas[i]) * this.omegas[i];
    vel.push([vx, vy]);
  }
  return vel;
};

Pendulum.prototype.getEnergy = function() {
  var energy = {kinetic: [], potential: []};
  var pos = this.getCartesianPositions();
  var vel = this.getCartesianVelocities();
  for (var i = 0; i < this.dimension; i++) {
    var t = 0.5 * this.masses[i] * (vel[i][0] * vel[i][0] + vel[i][1] * vel[i][1]);
    var u = this.masses[i] * GRAVITY * pos[i][1];
    energy.kinetic.push(t);
    energy.potential.push(u);
  }
  return energy;
};

Pendulum.prototype.getTotalEnergy = function() {
  var energy = this.getEnergy();
  var total_energy = 0;
  for (var i = 0; i < this.dimension; i++) {
    total_energy += energy.kinetic[i] + energy.potential[i];
  }
  return total_energy;
};

Pendulum.prototype.computeAlphas = function(thetas, omegas, torques) {
  if (this.dimension === 1) {
    return [-GRAVITY / this.lengths[0] * Math.sin(thetas[0])];
  } else if (this.dimension === 2) {
    var m1 = this.masses[0];
    var m2 = this.masses[1];
    var t1 = thetas[0];
    var t2 = thetas[1];
    var w1 = omegas[0];
    var w2 = omegas[1];
    var ell1 = this.lengths[0];
    var ell2 = this.lengths[1];
    
    var A = [[0, 0], [0, 0]];
    var b = [0, 0];

    A[0][0] = (m1 + m2) * ell1;
    A[0][1] = m2 * ell2 * Math.cos(t1 - t2);
    b[0] = -m2 * ell2 * w2 * w2 * Math.sin(t1 - t2);
    b[0] -= GRAVITY * (m1 + m2) * Math.sin(t1);
    A[1][0] = m2 * ell1 * Math.cos(t1 - t2);
    A[1][1] = m2 * ell2;
    b[1] = m2 * ell1 * ell2 * w1 * w1 * Math.sin(t1 - t2);
    b[1] -= ell2 * m2 * GRAVITY * Math.sin(t2);

    if (_.isArray(torques)) {
      b[0] += torques[0] * (t1 - t2) / ell1;
      b[1] -= torques[0] * (t1 - t2) / ell2;
    }

    return numeric.solve(A, b);
  } else {
    throw "Not implemented";
  }
}

// dt is in seconds
Pendulum.prototype.takeStep = function(dt, torques) {
  this.time += dt;
  if (this.integration_method === INTEGRATION_METHODS.FORWARD_EULER) {
    // console.log('forward')
    var alphas = this.computeAlphas(this.thetas, this.omegas);
    for (var i = 0; i < alphas.length; i++) {
      this.thetas[i] += this.omegas[i] * dt;
      this.omegas[i] += alphas[i] * dt;
    }
  } else if (this.integration_method === INTEGRATION_METHODS.SEMI_IMPLICIT_EULER) {
    // console.log('semi');
    var alphas = this.computeAlphas(this.thetas, this.omegas);
    for (var i = 0; i < alphas.length; i++) {
      this.omegas[i] += alphas[i] * dt;
      this.thetas[i] += this.omegas[i] * dt;
    }
  } else if (this.integration_method === INTEGRATION_METHODS.RUNGE_KUTTA) {
    var omegas_k1 = this.omegas.slice(0);
    var alphas_k1 = this.computeAlphas(this.thetas, this.omegas);

    var thetas = [];
    var omegas_k2 = [];
    for (var i = 0; i < alphas_k1.length; i++) {
      thetas.push(this.thetas[i] + (dt / 2) * omegas_k1[i]);
      omegas_k2.push(this.omegas[i] + (dt / 2) * alphas_k1[i]);
    }
    var alphas_k2 = this.computeAlphas(thetas, omegas_k2);

    var omegas_k3 = [];
    for (var i = 0; i < alphas_k2.length; i++) {
      thetas[i] = this.thetas[i] + (dt / 2) * omegas_k2[i];
      omegas_k3.push(this.omegas[i] + (dt / 2) * alphas_k2[i]);
    }
    var alphas_k3 = this.computeAlphas(thetas, omegas_k3);

    var omegas_k4 = [];
    for (var i = 0; i < alphas_k3.length; i++) {
      thetas[i] = this.thetas[i] + dt * omegas_k3[i];
      omegas_k4.push(this.omegas[i] + dt * alphas_k3[i]);
    }
    var alphas_k4 = this.computeAlphas(thetas, omegas_k4);

    for (var i = 0; i < this.dimension; i++) {
      var theta_step = omegas_k1[i] + 2 * omegas_k2[i] + 2 * omegas_k3[i] + omegas_k4[i];
      var omega_step = alphas_k1[i] + 2 * alphas_k2[i] + 2 * alphas_k3[i] + alphas_k4[i];
      this.thetas[i] += (dt / 6) * theta_step;
      this.omegas[i] += (dt / 6) * omega_step;
    }
  }
};

Pendulum.prototype.new_draw = function(canvas) {
  var ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  var scale = 0.5 * Math.min(canvas.width, canvas.height) / this.total_length;

  var pos = this.getCartesianPositions();
  for (var i = 0; i < this.dimension; i++) {
    if (i === 0) {
      var x0 = canvas.width / 2;
      var y0 = canvas.height / 2;
    } else {
      var x0 = canvas.width / 2 + scale * pos[i - 1][0];
      var y0 = canvas.height / 2 - scale * pos[i - 1][1];
    }
    var x1 = canvas.width / 2 + scale * pos[i][0];
    var y1 = canvas.height / 2 - scale * pos[i][1];
    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.stroke();

    ctx.beginPath();
    ctx.arc(x1, y1, scale * this.radii[i], 0, 2 * Math.PI);
    ctx.fill();
  }
};

// Draw this pendulum onto a canvas.
Pendulum.prototype.draw = function(canvas) {
  var ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // pixels / meter
  var scale = 0.5 * Math.min(canvas.width, canvas.height) / this.total_length;
  
  var x0 = canvas.width / 2;
  var y0 = canvas.height / 2;
  _.each(_.zip(this.lengths, this.masses, this.thetas, this.radii), function(elem, i) {
    var len = elem[0], mass = elem[1], theta = elem[2], radius = elem[3];

    var x1 = x0 + scale * Math.sin(theta) * len;
    var y1 = y0 + scale * Math.cos(theta) * len;

    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.stroke();

    var cx = x1 + scale * Math.sin(theta) * radius;
    var cy = y1 + scale * Math.cos(theta) * radius;
    ctx.beginPath();
    ctx.arc(x1, y1, scale * radius, 0, 2 * Math.PI);
    ctx.fill();

    x0 = x1;
    y0 = y1;
  });
};

var startSimulation = (function() {
  var intervals = [];
  return function(options) {
    _.each(intervals, clearInterval);
    intervals = [];
    var pendulum = new Pendulum(options);
    var canvas = $('#the-canvas').get(0);

    var energy_data = [];
    _.each(_.range(pendulum.dimension), function(i) {
      energy_data.push({label: 'Kinetic ' + (i + 1), data: []});
    });
    _.each(_.range(pendulum.dimension), function(i) {
      energy_data.push({label: 'Potential ' + (i + 1), data: []});
    });
    energy_data.push({label: 'Total energy', data: []});

    function update_energy_data() {
      var t = pendulum.time;
      var energy = pendulum.getEnergy();
      var total_energy = 0;
      for (var i = 0; i < pendulum.dimension; i++) {
        total_energy += energy.kinetic[i] + energy.potential[i];
        energy_data[i].data.push([t, energy.kinetic[i]]);
        energy_data[i + pendulum.dimension].data.push([t, energy.potential[i]]);
      }
      energy_data[2 * pendulum.dimension].data.push([t, total_energy]);
    }
    function draw_energy_plot() {
      var max_length = 1000;
      _.each(energy_data, function(series) {
        if (series.data.length > max_length) {
          series.data = series.data.slice(series.data.length - max_length);
        }
      });
      var flot_options = {
        legend: {position: 'nw'},
      };
      $.plot('#flot-div', energy_data, flot_options);
    }

    intervals.push(window.setInterval(function() {
      pendulum.takeStep(0.005);
      update_energy_data();
    }, 5));

    intervals.push(window.setInterval(function() {
      pendulum.new_draw(canvas);
    }, 10));

    draw_energy_plot();
    intervals.push(window.setInterval(draw_energy_plot, 100));
  }
})();

$(function() {
  $('#go-btn').click(function() {
    var options = {
      lengths: [1],
      masses: [2],
      thetas: [Math.PI / 2],
      omegas: [0],
      integration_method: INTEGRATION_METHODS.RUNGE_KUTTA,
    };
    var pendulum_type = $('#pendulum-type-btns input[type="radio"]:checked').val();

    if (pendulum_type === 'single') {
      // pass 
    } else if (pendulum_type === 'double') {
      options.lengths.push(1);
      options.masses.push(1);
      options.thetas.push(Math.PI / 2);
      options.omegas.push(0);
    }
    
    var integration_method = null;
    if ($('#int-forward').prop('checked')) {
      integration_method = INTEGRATION_METHODS.FORWARD_EULER;
    } else if ($('#int-semi').prop('checked')) {
      integration_method = INTEGRATION_METHODS.SEMI_IMPLICIT_EULER;
    } else if ($('#int-runge-kutta').prop('checked')) {
      integration_method = INTEGRATION_METHODS.RUNGE_KUTTA;
    }
    options.integration_method = integration_method;

    startSimulation(options);
  });
});
