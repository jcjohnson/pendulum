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
    'lengths', 'masses', 'thetas', 'omegas', 'integration_method', 'compound',
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
  _.each(_.zip(this.masses, this.lengths), function(x) {
    var m = x[0];
    var ell = x[1];
    if (this.compound) {
      this.radii.push(Math.sqrt(m / (Math.PI * ell * DENSITY)));
    } else {
      this.radii.push(Math.pow( (3 * m) / (4 * Math.PI * DENSITY), 1 / 3));
    }
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

Pendulum.prototype.getCentersOfMass = function() {
  if (!this.compound) return this.getCartesianPositions();
  var pos = [];
  var sx = 0;
  var sy = 0;
  for (var i = 0; i < this.dimension; i++) {
    var dsx = 0.5 * this.lengths[i] * Math.sin(this.thetas[i]);
    var dsy = 0.5 * this.lengths[i] * Math.cos(this.thetas[i]);
    sx += dsx;
    sy -= dsy;
    pos.push([sx, sy]);
    sx += dsx;
    sy -= dsy;
  }
  return pos;
};

Pendulum.prototype.getCartesianVelocities = function() {
  // TODO: update for compound
  var vel = [];
  var vx = 0;
  var vy = 0;
  for (var i = 0; i < this.dimension; i++) {
    var dvx = this.lengths[i] * Math.cos(this.thetas[i]) * this.omegas[i];
    var dvy = this.lengths[i] * Math.sin(this.thetas[i]) * this.omegas[i];
    if (this.compound) {
      dvx /= 2;
      dvy /= 2;
    }
    vx += dvx;
    vy += dvy;
    vel.push([vx, vy]);
    if (this.compound) {
      vx += dvx;
      vy += dvy;
    }
  }
  return vel;
};

Pendulum.prototype.getEnergyNames = function() {
  var range = _.range(this.dimension);
  var names = [];
  _.each(range, function(i) {
    names.push('Kinetic ' + (i + 1));
  });
  _.each(range, function(i) {
    names.push('Potential ' + (i + 1));
  });
  if (this.compound) {
    _.each(range, function(i) {
      names.push('Angular ' + (i + 1));
    });
  }
  names.push('Total energy');
  return names;
};

Pendulum.prototype.getEnergy = function() {
  var energy = [];
  var pos = this.getCentersOfMass();
  var vel = this.getCartesianVelocities();
  var total = 0;
  for (var i = 0; i < this.dimension; i++) {
    var t = 0.5 * this.masses[i] * (vel[i][0] * vel[i][0] + vel[i][1] * vel[i][1]);
    var u = this.masses[i] * GRAVITY * pos[i][1];
    total += t + u;
    energy.push(t);
    energy.push(u);
    if (this.compound) {
      var len = this.lengths[i];
      var w = this.omegas[i];
      var at = (1 / 24) * this.masses[i] * len * len * w * w;
      total += at;
      energy.push(at);
    }
  }
  energy.push(total);
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
    if (this.compound) {
      return [- (3 * GRAVITY) / (2 * this.lengths[0]) * Math.sin(thetas[0])];
    } else {
      return [-GRAVITY / this.lengths[0] * Math.sin(thetas[0])];
    }
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

    if (this.compound) {
      A[0][0] = (m1 / 3 + m2) * ell1;
      A[0][1] = 0.5 * m2 * ell2 * Math.cos(t1 - t2);
      b[0] = -0.5 * m2 * ell2 * w2 * w2 * Math.sin(t1 - t2);
      b[0] -= (m1 / 2 + m2) * GRAVITY * Math.sin(t1);
      A[1][0] = 0.5 * m2 * ell1 * Math.cos(t1 - t2);
      A[1][1] = (1 / 3) * m2 * ell2;
      b[1] = 0.5 * m2 * ell1 * w1 * w1 * Math.sin(t1 - t2);
      b[1] -= 0.5 * m2 * GRAVITY * Math.sin(t2);
    } else {
      A[0][0] = (m1 + m2) * ell1;
      A[0][1] = m2 * ell2 * Math.cos(t1 - t2);
      b[0] = -m2 * ell2 * w2 * w2 * Math.sin(t1 - t2);
      b[0] -= GRAVITY * (m1 + m2) * Math.sin(t1);
      A[1][0] = m2 * ell1 * Math.cos(t1 - t2);
      A[1][1] = m2 * ell2;
      b[1] = m2 * ell1 * ell2 * w1 * w1 * Math.sin(t1 - t2);
      b[1] -= ell2 * m2 * GRAVITY * Math.sin(t2);
    }

    /*
    // This adds a bit of stiffness
    b[0] -= 5 * (t1 - t2) / ell1;
    b[1] += 5 * (t1 - t2) / ell2;
    */

    /*
    if (_.isArray(torques)) {
      var tt1 = t1 % (2 * Math.PI);
      var tt2 = t2 % (2 * Math.PI);
      if (tt1 < 0) tt1 += 2 * Math.PI;
      if (tt2 < 0) tt2 += 2 * Math.PI;
      b[0] += (t1 - t2 - torques[0]) / ell1;
      b[1] -= (t1 - t2 - torques[0]) / ell2;
      // b[0] += torques[0] * (t1 - t2) / ell1;
      // b[1] -= torques[0] * (t1 - t2) / ell2;
    }
    */

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
    var alphas = this.computeAlphas(this.thetas, this.omegas, torques);
    for (var i = 0; i < alphas.length; i++) {
      this.thetas[i] += this.omegas[i] * dt;
      this.omegas[i] += alphas[i] * dt;
    }
  } else if (this.integration_method === INTEGRATION_METHODS.SEMI_IMPLICIT_EULER) {
    // console.log('semi');
    var alphas = this.computeAlphas(this.thetas, this.omegas, torques);
    for (var i = 0; i < alphas.length; i++) {
      this.omegas[i] += alphas[i] * dt;
      this.thetas[i] += this.omegas[i] * dt;
    }
  } else if (this.integration_method === INTEGRATION_METHODS.RUNGE_KUTTA) {
    var omegas_k1 = this.omegas.slice(0);
    var alphas_k1 = this.computeAlphas(this.thetas, this.omegas, torques);

    var thetas = [];
    var omegas_k2 = [];
    for (var i = 0; i < alphas_k1.length; i++) {
      thetas.push(this.thetas[i] + (dt / 2) * omegas_k1[i]);
      omegas_k2.push(this.omegas[i] + (dt / 2) * alphas_k1[i]);
    }
    var alphas_k2 = this.computeAlphas(thetas, omegas_k2, torques);

    var omegas_k3 = [];
    for (var i = 0; i < alphas_k2.length; i++) {
      thetas[i] = this.thetas[i] + (dt / 2) * omegas_k2[i];
      omegas_k3.push(this.omegas[i] + (dt / 2) * alphas_k2[i]);
    }
    var alphas_k3 = this.computeAlphas(thetas, omegas_k3, torques);

    var omegas_k4 = [];
    for (var i = 0; i < alphas_k3.length; i++) {
      thetas[i] = this.thetas[i] + dt * omegas_k3[i];
      omegas_k4.push(this.omegas[i] + dt * alphas_k3[i]);
    }
    var alphas_k4 = this.computeAlphas(thetas, omegas_k4, torques);

    for (var i = 0; i < this.dimension; i++) {
      var theta_step = omegas_k1[i] + 2 * omegas_k2[i] + 2 * omegas_k3[i] + omegas_k4[i];
      var omega_step = alphas_k1[i] + 2 * alphas_k2[i] + 2 * alphas_k3[i] + alphas_k4[i];
      this.thetas[i] += (dt / 6) * theta_step;
      this.omegas[i] += (dt / 6) * omega_step;
    }
  }
};

Pendulum.prototype.draw = function(canvas) {
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

    ctx.save();
    if (this.compound) {
      ctx.lineWidth = 2 * scale * this.radii[i];
    }

    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.stroke();
    ctx.restore();

    if (this.compound) {
      ctx.beginPath();
      ctx.arc(x0, y0, scale * this.radii[i], 0, 2 * Math.PI);
      ctx.fill();
    }

    ctx.beginPath();
    ctx.arc(x1, y1, scale * this.radii[i], 0, 2 * Math.PI);
    ctx.fill();
  }
};


var startSimulation = (function() {
  var intervals = [];
  return function(options) {
    _.each(intervals, clearInterval);
    intervals = [];
    var pendulum = new Pendulum(options);
    var canvas = $('#the-canvas').get(0);

    var energy_data = [];
    _.each(pendulum.getEnergyNames(), function(name) {
      energy_data.push({label: name, data: []});
    });

    function update_energy_data() {
      var t = pendulum.time;
      var energy = pendulum.getEnergy();
      for (var i = 0; i < energy.length; i++) {
        energy_data[i].data.push([t, energy[i]]);
      }
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
      var torque = [0];
      if (f_down) torque[0] += 0.5;
      if (j_down) torque[0] -= 0.5;
      pendulum.takeStep(0.005, torque);
      update_energy_data();
    }, 5));

    intervals.push(window.setInterval(function() {
      pendulum.draw(canvas);
    }, 10));

    draw_energy_plot();
    intervals.push(window.setInterval(draw_energy_plot, 100));

    var f_down = false;
    var j_down = false;

    $(document).off('keydown.pendulum');
    $(document).on('keydown.pendulum', function(e) {
      if (e.keyCode === 70) {
        f_down = true;
      } else if (e.keyCode === 74) {
        j_down = true;
      }
    });
    $(document).off('keyup.pendulum');
    $(document).on('keyup.pendulum', function(e) {
      if (e.keyCode === 70) {
        f_down = false;
      } else if (e.keyCode === 74) {
        j_down = false;
      }
      
    });
  }
})();

$(function() {
  $('#go-btn').click(function() {
    var options = {
      lengths: [1.2],
      masses: [40],
      thetas: [Math.PI / 2],
      omegas: [0],
      integration_method: INTEGRATION_METHODS.RUNGE_KUTTA,
      compound: true,
    };
    var pendulum_type = $('#pendulum-type-btns input[type="radio"]:checked').val();

    if (pendulum_type === 'single') {
      // pass 
    } else if (pendulum_type === 'double') {
      options.lengths.push(0.9);
      options.masses.push(40);
      options.thetas.push(Math.PI / 2);
      options.omegas.push(0);
    }

    if ($('#simple').prop('checked')) {
      options.compound = false;
    } else if ($('#compound').prop('checked')) {
      options.compound = true;
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
