// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// define constant parameters
// EPS is for the precision issue (see precept slide)
#define INFINITY 1.0e+12
#define EPS 1.0e-3

// define constants for scene setting 
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
    int shapeType;
    vec3 v1;
    vec3 v2;
    float rad;
};

struct Material {
    int materialType;
    vec3 color;
    float shininess;
    vec3 specular;

    int materialReflectType;
    float reflectivity; 
    float refractionRatio;
    int special;

};

struct Object {
    Shape shape;
    Material material;
};

struct Light {
    vec3 position;
    vec3 color;
    float intensity;
    float attenuate;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Intersection {
    vec3 position;
    vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;        
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

varying vec2 v_position;

// find the position some distance along a ray
vec3 rayGetOffset( Ray ray, float dist ) {
    return ray.origin + ( dist * ray.direction );
}

// if a newly found intersection is closer than the best found so far, record the new intersection and return true;
// otherwise leave the best as it was and return false.
bool chooseCloserIntersection( float dist, inout float best_dist, inout Intersection intersect, inout Intersection best_intersect ) {
    if ( best_dist <= dist ) return false;
    best_dist = dist;
    best_intersect.position = intersect.position;
    best_intersect.normal   = intersect.normal;
    return true;
}

// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 135 lines of code.
// float findTriangleArea(vec3 t1, vec3 t2, vec3 t3) {
//     vec3 e1 = t2 - t1;
//     vec3 e2 = t3 - t1;
//     vec3 e3 = cross(e1, e2);
//     float area = 0.5*length(e3);
//     return area;
// }

// float findNorm(vec3 t1, vec3 t2, vec3 t3, out vec3 norm) {
//     vec3 normal = cross(t2-t1, t3-t1);
//     norm = normalize(normal);
//     return dot(norm, t1 - t2);
// }

// bool checkPointInsideTriangle(vec3 p1, vec3 p2) {
// 	vec3 norm = normalize(cross( t2 - t1, t3 - t2 ));
//     float dist = dot( norm, t1 );
//     float len = findIntersectionWithPlane( ray, norm, dist, intersect );

//     vec3 v1 = t1 - intersect.position;
//     vec3 v2 = t2 - intersect.position;

//     vec3 n1 = normalize(cross( v1, v2 ));

//     if (dot( ray.direction, n1) < EPS) {
//     	return INFINITY;
//     }
// }

// float findDist(Ray ray, vec3 t1, vec3 t2, vec3 t3) {
    
// }
    // e1.x = t2.x - t1.x;
    // e1.y = t2.y - t1.y;
    // e1.z = t2.z - t1.z;

    // e2.x = t3.x - t1.x;
    // e2.y = t3.y - t1.y;
    // e2.z = t3.z - t3

// ----------- STUDENT CODE END ------------

// forward declaration
float rayIntersectScene( Ray ray, out Material out_mat, out Intersection out_intersect );

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane( Ray ray, vec3 norm, float dist, out Intersection intersect ) {
    float a   = dot( ray.direction, norm );
    float b   = dot( ray.origin, norm ) - dist;
    
    if ( a < EPS && a > EPS ) return INFINITY;
    
    float len = -b/a;
    if ( len < EPS ) return INFINITY;

    intersect.position = rayGetOffset( ray, len );
    intersect.normal   = norm;
    return len;
}
float findIntersectionFromThreePoints( Ray ray, vec3 t1, vec3 t2, vec3 t3, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 22 lines of code.
    
    vec3 norm = normalize(cross( t2 - t1, t3 - t2 ));
    float dist = dot( norm, t2 );
    return findIntersectionWithPlane( ray, norm, dist, intersect );
}
// Triangle
float findIntersectionWithTriangle( Ray ray, vec3 t1, vec3 t2, vec3 t3, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 22 lines of code.
    
    vec3 norm = normalize(cross( t2 - t1, t3 - t2 ));
    float dist = dot( norm, t2 );
    float len = findIntersectionWithPlane( ray, norm, dist, intersect );

    vec3 v1 = t1 - intersect.position;
    vec3 v2 = t2 - intersect.position;
    vec3 v3 = t3 - intersect.position;

    vec3 n1 = normalize(cross( v1, v2 ));

    if (dot( ray.direction, n1) < EPS) {
    	return INFINITY;
    }

    vec3 n2 = normalize(cross( v2, v3 ));

    if (dot( ray.direction, n2) < EPS) {
    	return INFINITY;
    }

    vec3 n3 = normalize(cross( v3, v1 ));

    if (dot( ray.direction, n3) < EPS) {
    	return INFINITY;
    }

    return len;
    
    // ----------- STUDENT CODE END ------------
}

// Sphere
float findIntersectionWithSphere( Ray ray, vec3 center, float radius, out Intersection intersect ) {   
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 23 lines of code.

    vec3 p0 = ray.origin;
    vec3 l = center - p0;
    vec3 v = ray.direction;
    float t_ca = dot(l, v);
    if (t_ca < EPS) {
    	return INFINITY;
    }
    float d_squared = dot(l, l) - (t_ca*t_ca);
    if (d_squared > radius*radius) {
    	return INFINITY;
    }
    float t_hc = sqrt(radius*radius - d_squared);
    float t_1 = t_ca - t_hc;
    float t_2 = t_ca + t_hc;

    if (t_1 > EPS) {
       	intersect.position = rayGetOffset(ray, t_1);
    	intersect.normal = normalize(intersect.position - center);
    	return t_1;
    }
    else if (t_2 > EPS) {
    	intersect.position = rayGetOffset(ray, t_2);
    	intersect.normal = normalize(intersect.position - center);
    	return t_2;
    }
    
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

// Box
float findIntersectionWithBox( Ray ray, vec3 pmin, vec3 pmax, out Intersection out_intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // pmin and pmax represent two bounding points of the box
    // pmin stores [xmin, ymin, zmin] and pmax stores [xmax, ymax, zmax]
    // ----------- Our reference solution uses 24 lines of code.
    // find three coordinates of front of box

    float t1;
    float t2;
    float t;
    float tmin = -INFINITY;
    float tmax = INFINITY;
    
    vec3 invRayDirection = 1.0 / ray.direction;
    if (ray.direction.x != EPS) {
        t1 = (pmin.x - ray.origin.x) * invRayDirection.x;
        t2 = (pmax.x - ray.origin.x) * invRayDirection.x;
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));
    }

    if (ray.direction.y != EPS) {
        t1 = (pmin.y - ray.origin.y) * invRayDirection.y;
        t2 = (pmax.y - ray.origin.y) * invRayDirection.y;
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));
    }

    if (ray.direction.z != EPS) {
        t1 = (pmin.z - ray.origin.z) * invRayDirection.z;
        t2 = (pmax.z - ray.origin.z) * invRayDirection.z;
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));
    } 

    if (tmin < EPS) {
        t = tmax;
    }
    else {
        t = tmin;
    }
    out_intersect.position = rayGetOffset( ray, t );
    return t;
    // vec3 t_bot = invRayDirection * (pmin - ray.origin);
    // vec3 t_top = invRayDirection * (pmax - ray.origin);
    // vec3 t_min = min(t_top, t_bot);
    // vec3 t_max = max(t_top, t_bot);
    // vec2 t = max(tmax.xx, tmin.yz);
    // float t_0 = max(t.x, t.y);
    // t = min(t_max.xx, t_max.yz);
    // float t_1 = min(t.x, t.y);
    // out_intersect.


    // vec3 s1 = vec3(pmin.x, pmin.y, pmin.z);
    // vec3 s2 = vec3(pmin.x, pmax.y, pmin.z);
    // vec3 s3 = vec3(pmax.x, pmin.y, pmin.z);

    // vec3 s4 = vec3(pmax.x, pmax.y, pmin.z);
    // vec3 s5 = vec3(pmax.x, pmax.y, pmax.z);
    // vec3 s6 = vec3(pmax.x, pmin.y, pmax.z);

    // vec3 s7 = vec3(pmin.x, pmin.y, pmax.z);
    // vec3 s8 = vec3(pmin.x, pmax.y, pmax.z);

    // float len1 = findIntersectionFromThreePoints( ray, s1, s2, s3, out_1);
    // float len2 = findIntersectionFromThreePoints( ray, s3, s4, s5, out_2 );

    // vec3 norm = normalize(cross( s2 - s1, s3 - s1 ));
    // float dist = dot( norm, s1 );
    // float len = findIntersectionWithPlane( ray, norm, dist, out_intersect );
    
    // // check if found intersection point is within the rectangle of the front of the box
    // if (out_intersect.position.x >= pmin.x && out_intersect.position.x <= pmax.x && out_intersect.position.y >= pmin.y && out_intersect.position.y <= pmax.y) {
    // 	return len;
    // }

    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}  

// Cylinder
float getIntersectOpenCylinder( Ray ray, vec3 center, vec3 axis, float len, float rad, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    // center = p_a + v_a*t
    // rearranging: p_a = center - v_a*t
    float t = INFINITY;
    vec3 v_a = axis;

    // vec3 p_0 = p_a;
    vec3 v = ray.direction;
    vec3 p = ray.origin;

    vec3 delta_p = p - center;

    vec3 a_1 = cross(v, v_a); //dot(cross(v, v_a), cross(v, v_a));    
    float a = dot(cross(v, v_a), cross(v, v_a));// dot((v - v_a*dot(v, v_a)),(v  - v_a*dot(v, v_a)));
    float b = 2.0 * dot(cross(v, v_a), cross(delta_p, v_a));//2*dot((v - dot(v, v_a)*v_a), delta_p - dot(delta_p, v_a)*v_a);
    float c = dot(cross(delta_p, v_a), cross(delta_p, v_a)) - (rad*rad * dot(v_a, v_a));//(delta_p - dot(delta_p, v_a)*v_a)*(delta_p - dot(delta_p, v_a)*v_a) - r*r;
    // float dotDirAxis = dot(v, v_a);
    // vec3 aVec = v - dotDirAxis*v_a;

    // float a = dot(aVec, aVec);
    // float b = 2.0*dot(aVec, delta_p - dot(delta_p, v_a)*v_a);
    // float c = dot(delta_p - dot(delta_p, v_a)*v_a, delta_p - dot(delta_p, v_a)*v_a) - rad*rad;

    float b24ac = sqrt(b*b - 4.0*a*c);
    // if (b24ac < EPS) return INFINITY;
    float tInt1 = (-b + b24ac)/(2.0*a);
    float tInt2 = (-b - b24ac)/(2.0*a);

    vec3 center2 = center + axis * len;
    // float tInt1 = (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    // float tInt2 = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

    float tInt = tInt1;
    vec3 intPt = rayGetOffset(ray, tInt);
    float heightCoord = dot(axis, intPt - center);
    if (tInt < EPS || heightCoord < 0.0 || heightCoord > len) {
        tInt = tInt2;
        intPt = rayGetOffset(ray, tInt);
        heightCoord = dot(axis, intPt - center);
        if (tInt < EPS || heightCoord < 0.0 || heightCoord > len) {
            return INFINITY;
        }
    }

    vec3 cylinderCenterAtHeight = center + axis * heightCoord;
    intersect.position = intPt;
    intersect.normal = normalize(intPt - cylinderCenterAtHeight);
    return tInt;

    // vec3 q_0 = rayGetOffset(ray, t_0);
    // vec3 q_1 = rayGetOffset(ray, t_1);


    // t = t_0;
    // vec3 q = q_0;
    // float heightCoord = dot(v_a, q - center);

    // if (t < EPS || heightCoord < EPS || heightCoord > len) {
    //     t = t_1;
    //     q = q_1;
    //     heightCoord = dot(v_a, q - center);
    //     if (t < EPS || heightCoord < 0.0 || heightCoord > len) {
    //         return INFINITY;
    //     }
    // }

    // // if (t_0 > EPS || t_1 > EPS) {
    // //     // dot(v_a, (q_0 - p_0)) is the projection of v_a along the  
    // //     if (t_0 > EPS && t_0 < ) {//&& dot(v_a, (q_0 - p_0)) > EPS && dot(v_a, (q_0 - p_1)) < EPS) {
    // //         t = t_0;
    // //     }
    // //     if (t_1 > EPS && t_1 < 1.0) {//dot(v_a, (q_1 - p_0)) > EPS && dot(v_a, (q_1 - p_1)) < EPS) {
    // //         if (t_1 < t_0) {
    // //             t = t_1;

    // //         }
    // //     }

    
    // // vec3 cylinderCenterAtHeight = center + axis * heightCoord;

    // // intersect.position = rayGetOffset( ray, t );
    // // intersect.normal = normalize(intersect.position - cylinderCenterAtHeight);
    // // return t;

    // vec3 cylinderCenterAtHeight = center + axis * heightCoord;
    // intersect.position = q;
    // intersect.normal = normalize(t - cylinderCenterAtHeight);
    // return t;

    


 // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

float getIntersectDisc( Ray ray, vec3 center, vec3 norm, float rad, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 15 lines of code
    float dist = dot(norm, center);

    float len = findIntersectionWithPlane( ray, norm, dist, intersect );
    vec3 v = center - intersect.position;
    float d_squared = dot(v, v);

    if ( d_squared <= rad*rad) {
    	return len;
    }
    else {
    	return INFINITY;
    }
    // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}


float findIntersectionWithCylinder( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis = apex - center;
    float len = length( axis );
    axis = normalize( axis );

    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cylinder
    dist = getIntersectOpenCylinder( ray, center, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- two caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );
    dist = getIntersectDisc( ray,   apex, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // dist = findIntersectionWithPlane( ray, axis, radius, intersect );
    // chooseCloserIntersection( dist, best_dist, intersect, out_intersect );
    // dist = findIntersectionWithPlane( ray, axis, radius, intersect );
    // chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}
    
// Cone
float getIntersectOpenCone( Ray ray, vec3 apex, vec3 axis, float len, float radius, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    return INFINITY; // currently reports no intersection
    // ----------- STUDENT CODE END ------------
}

float findIntersectionWithCone( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis   = center - apex;
    float len   = length( axis );
    axis = normalize( axis );
        
    // -- infinite cone
    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cone
    dist = getIntersectOpenCone( ray, apex, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}

#define MAX_RECURSION 8

vec3 calculateSpecialDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    // ----------- STUDENT CODE BEGIN ------------
    if ( mat.special == CHECKERBOARD ) {
        // do something here for checkerboard
        // ----------- Our reference solution uses 21 lines of code.
        if (mod(floor(posIntersection.x) + floor(posIntersection.y + floor(posIntersection.z)), float(2)) == 1.0) {
            return vec3(1.0, 1.0, 1.0);
            
        }
        else {
            return vec3(EPS, EPS, EPS);
        }
    } 
    else if ( mat.special == MYSPECIAL ) {
        // do something here for myspecial
        // ----------- Our reference solution uses 2 lines of code.
    }

    return mat.color; // special materials not implemented. just return material color.
    // ----------- STUDENT CODE END ------------
}

vec3 calculateDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    // Special colors
    if ( mat.special != NONE ) {
        return calculateSpecialDiffuseColor( mat, posIntersection, normalVector ); 
    }
    return vec3( mat.color );
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light
bool pointInShadow( vec3 pos, vec3 lightVec ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 10 lines of code.
    return false;
    // ----------- STUDENT CODE END ------------
}

vec3 getLightContribution( Light light, Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly, vec3 diffuseColor ) {

    vec3 lightVector = light.position - posIntersection;
    
    if ( pointInShadow( posIntersection, lightVector ) ) {
        return vec3( EPS, EPS, EPS );
    }

    if ( mat.materialType == PHONGMATERIAL || mat.materialType == LAMBERTMATERIAL ) {
        vec3 contribution = vec3( EPS, EPS, EPS );

        // get light attenuation
        float dist = length( lightVector );
        float attenuation = light.attenuate * dist * dist;

        float diffuseIntensity = max( EPS, dot( normalVector, lightVector ) ) * light.intensity;
        
        // glass and mirror objects have specular highlights but no diffuse lighting
        if ( !phongOnly ) {
            contribution += diffuseColor * diffuseIntensity * light.color / attenuation;
        }
        
        if ( mat.materialType == PHONGMATERIAL ) {
            // ----------- STUDENT CODE BEGIN ------------
            vec3 phongTerm = vec3( EPS, EPS, EPS ); // not implemented yet, so just add black   
            // ----------- Our reference solution uses 10 lines of code.
            // ----------- STUDENT CODE END ------------
            contribution += phongTerm;
        }

        return contribution;
    }
    else {
        return diffuseColor;
    }

}

vec3 calculateColor( Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly ) {
	vec3 diffuseColor = calculateDiffuseColor( mat, posIntersection, normalVector );

	vec3 outputColor = vec3( EPS, EPS, EPS ); // color defaults to black when there are no lights
	
    for ( int i=0; i<MAX_LIGHTS; i++ ) {

        if( i>=numLights ) break; // because GLSL will not allow looping to numLights
		
        outputColor += getLightContribution( lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor );
	}
	
	return outputColor;
}

// find reflection or refraction direction ( depending on material type )
vec3 calcReflectionVector( Material material, vec3 direction, vec3 normalVector, bool isInsideObj ) {
    if( material.materialReflectType == MIRRORREFLECT ) {
        return reflect( direction, normalVector );
    }
    // the material is not mirror, so it's glass.
    // compute the refraction direction...
    
    // ----------- STUDENT CODE BEGIN ------------
    // see lecture 13 slide ( lighting ) on Snell's law
    // the eta below is eta_i/eta_r
    float eta = ( isInsideObj ) ? 1.0/material.refractionRatio : material.refractionRatio;
    // ----------- Our reference solution uses 11 lines of code.
    
    return reflect( direction, normalVector ); // return mirror direction so you can see something
    // ----------- STUDENT CODE END ------------
}

vec3 traceRay( Ray ray ) {
    Material hitMaterial;
    Intersection intersect;

    vec3 resColor  = vec3( EPS, EPS, EPS );
    vec3 resWeight = vec3( 1.0, 1.0, 1.0 );
    
    bool isInsideObj = false;

    for ( int depth = 0; depth < MAX_RECURSION; depth++ ) {
        
        float hit_length = rayIntersectScene( ray, hitMaterial, intersect );
            
        if ( hit_length < EPS || hit_length >= INFINITY ) break;

        vec3 posIntersection = intersect.position;
        vec3 normalVector    = intersect.normal;

        vec3 eyeVector = normalize( ray.origin - posIntersection );           
        if ( dot( eyeVector, normalVector ) < EPS )
            { normalVector = -normalVector; isInsideObj = true; }
        else isInsideObj = false;

        bool reflective = ( hitMaterial.materialReflectType == MIRRORREFLECT || 
                            hitMaterial.materialReflectType == GLASSREFLECT );
		vec3 outputColor = calculateColor( hitMaterial, posIntersection, normalVector, eyeVector, reflective );

        float reflectivity = hitMaterial.reflectivity;

        // check to see if material is reflective ( or refractive )
        if ( !reflective || reflectivity < EPS ) {
            resColor += resWeight * outputColor;
            break;
        }
        
        // bounce the ray
        vec3 reflectionVector = calcReflectionVector( hitMaterial, ray.direction, normalVector, isInsideObj );
        ray.origin = posIntersection;
        ray.direction = normalize( reflectionVector );

        // add in the color of the bounced ray
        resColor += resWeight * outputColor;
        resWeight *= reflectivity;
    }

    return resColor;
}

void main( ) {
    float cameraFOV = 0.8;
    vec3 direction = vec3( v_position.x * cameraFOV * width/height, v_position.y * cameraFOV, 1.0 );

    Ray ray;
	ray.origin    = vec3( uMVMatrix * vec4( camera, 1.0 ) );
    ray.direction = normalize( vec3( uMVMatrix * vec4( direction, EPS ) ) );

    // trace the ray for this pixel
    vec3 res = traceRay( ray );
    
    // paint the resulting color into this pixel
    gl_FragColor = vec4( res.x, res.y, res.z, 1.0 );
}

