function matxvec(matrix, vector)
{
    let result = []
    if (matrix.length == 3)
    {
        result.push(matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2]);
        result.push(matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2]);
        result.push(matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2]);
    }

    if (matrix.length == 9)
    {
        result.push(matrix[0] * vector[0] + matrix[1] * vector[1] + matrix[2] * vector[2]);
        result.push(matrix[3] * vector[0] + matrix[4] * vector[1] + matrix[5] * vector[2]);
        result.push(matrix[6] * vector[0] + matrix[7] * vector[1] + matrix[8] * vector[2]);
    }

    return result;
}

function sumvecs(a, b)
{
    let result = [];
    a.forEach(function (e, idx) {
        result.push(e + b[idx]);
    });
    return result;
}

function sumarrayofvecs(array)
{
    let result = [0, 0, 0];
    for (let vec = 0; vec < array.length; vec++)
    {
        result = sumvecs(result, array[vec]);
    }

    return result;
}

function diag(a)
{
    return [a, 0, 0, 0, a, 0, 0, 0, a];
}

function difvecs(a, b)
{
    let result = [];
    a.forEach(function (e, idx) {
        result.push(e - b[idx]);
    });
    return result;
}

function bdmxvec(bdm, vector)
{
    let result = [];

    for (let i = 0; i < bdm.length; i++)
        result.push(matxvec(bdm[i], vector[i]));
    
    return result;
}

function numxvec(num, vector)
{
    let result = [];
    vector.forEach(function (e) {
        result.push(num * e);
    });
    return result;
}

function dot(a, b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function length(x)
{
    return Math.sqrt(dot(x, x));
}

function cross(a, b)
{
    return ([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]);
}

//--------------------------------------------------taken from MDN

multiplyMatrices = function (a, b) {
  
    
    var result = [];
    
    var a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
        a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
        a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
        a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];

    var b0  = b[0], b1 = b[1], b2 = b[2], b3 = b[3];  
    result[0] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
    result[1] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
    result[2] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
    result[3] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
  
    b0 = b[4]; b1 = b[5]; b2 = b[6]; b3 = b[7];
    result[4] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
    result[5] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
    result[6] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
    result[7] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
  
    b0 = b[8]; b1 = b[9]; b2 = b[10]; b3 = b[11];
    result[8] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
    result[9] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
    result[10] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
    result[11] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
  
    b0 = b[12]; b1 = b[13]; b2 = b[14]; b3 = b[15];
    result[12] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
    result[13] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
    result[14] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
    result[15] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
  
    return result;
  }
  
  multiplyArrayOfMatrices = function (matrices) {
    
    var inputMatrix = matrices[0];
    
    for(var i=1; i < matrices.length; i++) {
      inputMatrix = multiplyMatrices(inputMatrix, matrices[i]);
    }
    
    return inputMatrix;
  }
  
  rotateXMatrix = function (a) {
    
    var cos = Math.cos;
    var sin = Math.sin;
    
    return [
         1,       0,        0,     0,
         0,  cos(a),  -sin(a),     0,
         0,  sin(a),   cos(a),     0,
         0,       0,        0,     1
    ];
  }
  
  rotateYMatrix = function (a) {
  
    var cos = Math.cos;
    var sin = Math.sin;
    
    return [
       cos(a),   0, sin(a),   0,
            0,   1,      0,   0,
      -sin(a),   0, cos(a),   0,
            0,   0,      0,   1
    ];
  }
  
  rotateZMatrix = function (a) {
  
    var cos = Math.cos;
    var sin = Math.sin;
    
    return [
      cos(a), -sin(a),    0,    0,
      sin(a),  cos(a),    0,    0,
           0,       0,    1,    0,
           0,       0,    0,    1
    ];
  }
  
  translateMatrix = function (x, y, z) {
      return [
          1,    0,    0,   0,
          0,    1,    0,   0,
          0,    0,    1,   0,
          x,    y,    z,   1
      ];
  }
  
  scaleMatrix = function (w, h, d) {
      return [
          w,    0,    0,   0,
          0,    h,    0,   0,
          0,    0,    d,   0,
          0,    0,    0,   1
      ];
  }
  
  perspectiveMatrix = function (fieldOfViewInRadians, aspectRatio, near, far) {
    
    // Construct a perspective matrix
    
    /*
       Field of view - the angle in radians of what's in view along the Y axis
       Aspect Ratio - the ratio of the canvas, typically canvas.width / canvas.height
       Near - Anything before this point in the Z direction gets clipped (outside of the clip space)
       Far - Anything after this point in the Z direction gets clipped (outside of the clip space)
    */
    
    var f = 1.0 / Math.tan(fieldOfViewInRadians / 2);
    var rangeInv = 1 / (near - far);
   
    return [
      f / aspectRatio, 0,                          0,   0,
      0,               f,                          0,   0,
      0,               0,    (near + far) * rangeInv,  -1,
      0,               0,  near * far * rangeInv * 2,   0
    ];
  }
  
  orthographicMatrix = function(left, right, bottom, top, near, far) {
    
    // Each of the parameters represents the plane of the bounding box
    
    var lr = 1 / (left - right);
    var bt = 1 / (bottom - top);
    var nf = 1 / (near - far);
      
    var row4col1 = (left + right) * lr;
    var row4col2 = (top + bottom) * bt;
    var row4col3 = (far + near) * nf;
    
    return [
       -2 * lr,        0,        0, 0,
             0,  -2 * bt,        0, 0,
             0,        0,   2 * nf, 0,
      row4col1, row4col2, row4col3, 1
    ];
  }