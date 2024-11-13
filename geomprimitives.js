//Purpose: The engine behind the 3D primitive operations for Mini Assignment 1

vec3 = glMatrix.vec3;

//////////////////////////////////////////////
///********         PART 1          *******///
//////////////////////////////////////////////


/**
 * 
 * @param {vec3} u
 * @returns {float} The magnitude of vector u
 */
function getMagnitude(u){
    return Math.sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])
}

/**
 * Compute the angle between the vectors ab and ac
 * @param {vec3} a First point
 * @param {vec3} b Second point
 * @param {vec3} c Third point
 * 
 * @return {float} Angle between vectors ab and ac in degrees
 */
function getAngle(a, b, c) {
    // TODO: Fill this in
    let ret = -1

    let ab = vec3.subtract(vec3.create(), b, a)
    let ac = vec3.subtract(vec3.create(), c, a)

    if (getMagnitude(ab) > 0 && getMagnitude(ac) > 0){
        ret = Math.acos(vec3.dot(ab, ac) / (getMagnitude(ab) * getMagnitude(ac))) * 180/Math.PI
    }

    return ret
}



/**
 * Project vector u onto vector v using the glMatrix library
 * @param {vec3} u Vector that's being projected
 * @param {vec3} v Vector onto which u is projected
 * 
 * @return {vec3} The projection of u onto v
 */
function projVector(u, v) {
    // TODO: Fill this in
    let ret = vec3.fromValues(0,0,0)
    if (getMagnitude(v) > 0){
        ret = vec3.scale(vec3.create(), v, vec3.dot(u, v) / vec3.dot(v, v))
    }
    return ret
}


/**
 * 
 * @param {vec3} u Vector that's being projected
 * @param {vec3} v Vector onto which u is perpendicularly projected
 * 
 * @return {vec3} The perpendicular projection of u onto v
 */
function projPerpVector(u, v) {
    // TODO: Fill this in
    let ret = vec3.fromValues(0,0,0)
    if (getMagnitude(v) > 0){
        ret = vec3.subtract(vec3.create(), u, projVector(u, v))
    }
    return ret
}


/**
 * Given three 3D vertices a, b, and c, compute the area 
 * of the triangle they span
 * @param {vec3} a First point
 * @param {vec3} b Second point
 * @param {vec3} c Third point
 * 
 * @return {float} Area of the triangle
 */
function getTriangleArea(a, b, c) {
    // TODO: Fill this in
    let ab = vec3.subtract(vec3.create(), b, a)
    let ac = vec3.subtract(vec3.create(), c, a)
    let crossProduct = vec3.cross(vec3.create(), ab, ac)
    return getMagnitude(crossProduct)/2
}


/**
 * For a plane determined by the points a, b, and c, with the plane
 * normal determined by those points in counter-clockwise order using 
 * the right hand rule, decide whether the point d is above, below, or on the plane
 * @param {vec3} a First point on plane
 * @param {vec3} b Second point on plane
 * @param {vec3} c Third point on plane
 * @param {vec} d Test point
 * 
 * @return {int} 1 if d is above, -1 if d is below, 0 if d is on
 */
function getAboveOrBelow(a, b, c, d) {
    // TODO: Fill this in
    let ab = vec3.subtract(vec3.create(), b, a)
    let ac = vec3.subtract(vec3.create(), c, a)
    let normal = vec3.cross(vec3.create(), ab, ac)
    let ad = vec3.subtract(vec3.create(), d, a)
    let dotproduct = vec3.dot(ad, normal)
    let ret

    const isCollinear = (a, b, c) => {
        let ab = vec3.subtract(vec3.create(), b, a)
        let ac = vec3.subtract(vec3.create(), c, a)
        let bc = vec3.subtract(vec3.create(), c, b)

        let one = getMagnitude(vec3.cross(vec3.create(), ab, ac))
        let two = getMagnitude(vec3.cross(vec3.create(), ab, bc))
        let three = getMagnitude(vec3.cross(vec3.create(), ac, bc))
        
        return one === 0 && two === 0 && three === 0
    }

    if (isCollinear(a, b, c)){
        ret = -2
    } else if (dotproduct < 0){
        ret = -1
    } else if (dotproduct > 0){
        ret = 1
    } else {
        ret = 0
    }
    return ret
}


//////////////////////////////////////////////
///********         PART 2          *******///
//////////////////////////////////////////////

/**
 * Compute the barycentric coordinates of a point p with respect to a triangle /\abc
 * 
 * @param {vec3} a Point a on the triangle
 * @param {vec3} b Point b on the triangle
 * @param {vec3} c Point c on the triangle
 * @param {vec3} p The point whose barycentric coordinates we seek
 * 
 * @return {vec3} An vec3 with the barycentric coordinates (alpha, beta, gamma)
 * 				  corresponding to a, b, and c, respectively, so that
 * 				  alpha + beta + gamma = 1, and alpha, beta, gamma >= 0
 *          CORNER CASES:
 * 				  (1) If p is not inside of /\abc, then return [0, 0, 0]
 *          (2) If /\abc is zero area, then return [1, 0, 0] iff p = a (=b=c)
 *              otherwise, return [0, 0, 0]
 */
function getBarycentricCoords(a, b, c, p) {
    // TODO: Fill this in
    let ret
    let A = getTriangleArea(a, b, c)
    let Aa = getTriangleArea(b, c, p)
    let Ab = getTriangleArea(a, c, p)
    let Ac = getTriangleArea(a, b, p)

    let alpha = Aa/A
    let beta = Ab/A
    let gamma = Ac/A

    if (A === 0 && vec3.equals(a, p)){
        ret = [1,0,0]
    } else if (A === 0 || alpha + beta + gamma > 1 + 1e-5){
        ret = [0,0,0]
    }else {
        ret = vec3.fromValues(alpha, beta, gamma)
    }

    return ret
}


/**
 * Find the intersection of a ray with a triangle
 * 
 * @param {vec3} p0 Endpoint of ray 
 * @param {vec3} v Direction of ray
 * @param {vec3} a Triangle vertex 1
 * @param {vec3} b Triangle vertex 2
 * @param {vec3} c Triangle vertex 3
 * 
 * @return {list} A list of vec3 objects.  The list should be empty
 *          if there is no intersection, or it should contain 
 *          exactly one vec3 object if there is an intersection
 *          CORNER CASES:
 *          (1) If the ray is parallel to the plane,
*               return an empty list
 */
function rayIntersectTriangle(p0, v, a, b, c) {
	// TODO: Fill this in
    let ret = []

    let ab = vec3.subtract(vec3.create(), b, a)
    let ac = vec3.subtract(vec3.create(), c, a)
    let n = vec3.cross(vec3.create(), ab, ac)

    // plane intersection
    let t = ( vec3.dot( vec3.subtract(vec3.create(), a, p0), n ) / vec3.dot(v, n) )
    let plane_intersection = vec3.scaleAndAdd(vec3.create(), p0, v, t)

    let coords = getBarycentricCoords(a, b, c, plane_intersection)

    inTriangle = coords[0]+coords[1]+coords[2] < 1+1e-5 && coords[0]+coords[1]+coords[2] > 1-1e-5

    if (inTriangle && vec3.dot(v, n) != 0 && t >= 0){
        ret = [plane_intersection]
    }


	return ret
}


/**
 * Find the intersection of the ray p0 + tv, t >= 0, with the
 * sphere centered at c with radius r.
 * 
 * @param {vec3} p0 Endpoint of the ray
 * @param {vec3} v Direction of the ray
 * @param {vec3} c Center of the sphere
 * @param {number} r Radius of the sphere
 * 
 * @return {list of vec3} A list of intersection points, 
 *   ***in the order in which the ray hits them***
 * If the ray doesn't hit any points, this list should be empty.
 * Note that a ray can hit at most 2 points on a sphere.
 */
function rayIntersectSphere(p0, v, center, radius) {
	// TODO: Fill this in
    let ret
    let w = vec3.subtract(vec3.create(), p0, center)
    v = vec3.normalize(v, v)

    let a = vec3.dot(v, v)
    let b = vec3.dot(vec3.scale(vec3.create(), w, 2), v)
    let c = vec3.dot(w, w) - radius*radius

    let discriminant = b*b - 4*a*c
    
    if (discriminant < 0){
        ret = []
    } else if (discriminant === 0){
        t = -1*b/2*a
        ret = [vec3.scaleAndAdd(vec3.create(), p0, v, t)]
    } else {
        // find solutions and check for special cases
        let tmin = (-1*b - Math.sqrt(discriminant))/2*a
        let tmax = (-1*b + Math.sqrt(discriminant))/2*a

        if (tmin < 0 && tmax < 0){
            ret = []
        } else if (tmin < 0 && tmax >= 0){
            ret = [vec3.scaleAndAdd(vec3.create(), p0, v, tmax)]
        } else if (tmin >= 0 && tmax < 0){
            ret = [vec3.scaleAndAdd(vec3.create(), p0, v, tmin)]
        } else {
            ret = [vec3.scaleAndAdd(vec3.create(), p0, v, tmin), vec3.scaleAndAdd(vec3.create(), p0, v, tmax)]
        }
    }
    return ret
}