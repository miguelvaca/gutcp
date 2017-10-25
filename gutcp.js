/*
 * File   : gutcp.js
 * Author : J Miguel Vaca
 * Date   : 04-OCT-2017
 * Contains the reusable classes for the GUTCP-based visualisations/simulations. Requires THREE.js to be loaded 
 * beforehand.
 */

var cvfVertexShader = 
	"uniform float wavesPerRing; 	\n" + 
	"uniform float waveSpeed; 		\n" + 
	"uniform float uTime; 			\n" + 
	"uniform float scale; 			\n" + 
	"uniform bool  animate; 		\n" + 
	"uniform bool  swapAnimDir; 	\n" + 
	"uniform vec4  color; 			\n" + 
	"uniform vec4  background; 		\n" + 
	"varying vec4  vColor; 			\n" + 
	"\n" 								+
	"void main() { 					\n" + 
	"	if(animate) { 				\n" + 
	"		float bright = mod(uTime + (swapAnimDir ? uv.y : uv.x), 1.0/wavesPerRing) * wavesPerRing; \n" + 
	"		//vColor = vec4( (1.0-bright)*color[0], (1.0-bright)*color[1], (1.0-bright)*color[2], color[3] ); \n" + 
	"		vColor = vec4( (1.0-bright)*color[0]+bright*background[0], (1.0-bright)*color[1]+bright*background[1], (1.0-bright)*color[2]+bright*background[2], color[3] );	\n" + 
	"	} else { 																					\n" + 
	"		vColor = color; 																		\n" + 
	"	} 																							\n" + 
	"	\n" 																							+	
	"	vec4 mvPosition = modelViewMatrix * vec4( position.xyz*scale, 1.0 ); 						\n" + 
	"	gl_PointSize = 5.0; 																		\n" + 
	"	gl_Position = projectionMatrix * mvPosition; 												\n" + 
	"} \n";

var cvfVertexShader2 = 
	"uniform float wavesPerRing; 	\n" + 
	"uniform float waveSpeed; 		\n" + 
	"uniform float uTime; 			\n" + 
	"uniform float scale; 			\n" + 
	"uniform float sphere; 			\n" + 
	"uniform bool  animate; 		\n" + 
	"uniform bool  swapAnimDir; 	\n" + 
	"uniform vec3  pole; 			\n" + 
	"uniform vec4  color; 			\n" + 
	"varying vec4  vColor; 			\n" + 
	"const float PI = 3.1415926535897932384626433832795; \n" + 
	"\n" 								+
	"void main() { 					\n" + 
	"   //sphere = (sphere < 0.01) ? 0.01 : sphere; \n" + 
	"   vec3  sphere_center = vec3(0.0, 0.0, 0.0); \n" + 
	"   //vec3  pole    = vec3(100.0*sin(pole_theta*0.01745)*sin(pole_psi*0.01745), 100.0*cos(pole_theta*0.01745), 100.0*sin(pole_theta*0.01745)*cos(pole_psi*0.01745)); \n" + 
	"   //float radius  = distance(sphere_center, pole); \n" + 
	"   float radius  = length(pole); \n" + 
	"   vec3  sc_pole = pole - sphere_center; \n" + 
	"   vec3  sc_pos  = position - sphere_center; \n" + 
	"   float initial_angle = acos( dot(sc_pole, sc_pos) / (length(sc_pole) * length(sc_pos)) ); \n" + 
	"   float arc_length = initial_angle * radius; \n" + 
	"   vec3  new_sphere_center = (-sc_pole / sphere) + pole; \n" + 
	"   float new_radius = radius / sphere; \n" + 
	"   float new_angle  = initial_angle * sphere; \n" + 
	"   vec3  new_sc_pole = pole - new_sphere_center; \n" + 
	"   vec3  new_sc_pos  = position - new_sphere_center; \n" + 
	"   vec3  rotation_axis = normalize( cross(new_sc_pole, new_sc_pos) ); \n" + 
	"   float cos_new_angle = cos(new_angle); \n" + 
	"   float sin_new_angle = sin(new_angle); \n" + 
	"   mat3  rotation_m = mat3( cos_new_angle+(1.0-cos_new_angle)*pow(rotation_axis.x,2.0), 						\n" + 
	"                            (1.0-cos_new_angle)*rotation_axis.x*rotation_axis.y+sin_new_angle*rotation_axis.z, \n" + 
	"                            (1.0-cos_new_angle)*rotation_axis.x*rotation_axis.z-sin_new_angle*rotation_axis.y, \n" +
	"                            (1.0-cos_new_angle)*rotation_axis.x*rotation_axis.y-sin_new_angle*rotation_axis.z, \n" +
	"                            cos_new_angle+(1.0-cos_new_angle)*pow(rotation_axis.y,2.0), 						\n" + 
	"                            (1.0-cos_new_angle)*rotation_axis.y*rotation_axis.z+sin_new_angle*rotation_axis.x, \n" +
	"                            (1.0-cos_new_angle)*rotation_axis.x*rotation_axis.z+sin_new_angle*rotation_axis.y, \n" +
	"                            (1.0-cos_new_angle)*rotation_axis.y*rotation_axis.z-sin_new_angle*rotation_axis.x, \n" +
	"                            cos_new_angle+(1.0-cos_new_angle)*pow(rotation_axis.z,2.0) );  					\n" + 
	"   vec3  new_position = rotation_m * new_sc_pole + sphere_center; \n" + 
	"   //vec3  new_position = rotation_m * (length(new_sc_pos) * normalize(new_sc_pole)) + new_sphere_center; \n" + 
	"   \n" + 
	"   \n" + 
	"	if(animate) { 				\n" + 
	"	float bright = mod(waveSpeed*uTime + (swapAnimDir ? uv.y : uv.x), 1.0/wavesPerRing) * wavesPerRing; \n" + 
	"		vColor = vec4( (1.0-bright)*color[0], (1.0-bright)*color[1], (1.0-bright)*color[2], color[3] );	\n" + 
	"	} else { 																							\n" + 
	"		vColor = color; 																				\n" + 
	"	} 																									\n" + 
	"	\n" 																									+	
	"	vec4 mvPosition = modelViewMatrix * vec4( new_position.xyz*scale, 1.0 ); 							\n" + 
	"	gl_PointSize = 5.0; 																				\n" + 
	"	gl_Position = projectionMatrix * mvPosition; 														\n" + 
	"} \n";
/*
		if(mode == 1) {
			// Now rotate RH about (-1, 1, 0) axis by -90 degrees, so that it can move as a TEM wave in the +z direction:
			// [Or equivalently, about (1, -1, 0) by 90 degrees.]
			vertex2[ 0 ] =    0.5*vertex[0] -   0.5*vertex[1] - 0.707*vertex[2];
			vertex2[ 1 ] =   -0.5*vertex[0] +   0.5*vertex[1] - 0.707*vertex[2]; 
			vertex2[ 2 ] =  0.707*vertex[0] + 0.707*vertex[1] + 0.0;
		} else {
			// Now rotate LH about (1, 1, 0) axis by +90 degrees, so that it can move as a TEM wave in the +z direction:
			vertex2[ 0 ] =    0.5*vertex[0] +   0.5*vertex[1] + 0.707*vertex[2];
			vertex2[ 1 ] =    0.5*vertex[0] +   0.5*vertex[1] - 0.707*vertex[2]; 
			vertex2[ 2 ] = -0.707*vertex[0] + 0.707*vertex[1] + 0.0;
		}
*/
var photonVertexShader = 
	"uniform float wavesPerRing; 	\n" + 
	"uniform float waveSpeed; 		\n" + 
	"uniform float uTime; 			\n" + 
	"uniform float radius; 			\n" + 
	"uniform float scale; 			\n" + 
	"uniform float cutoff; 			\n" + 
	"uniform bool  animate; 		\n" + 
	"uniform bool  swapAnimDir; 	\n" + 
	"uniform bool  rotate; 			\n" + 
	"uniform bool  rhcp; 			\n" + 
	"uniform vec4  color; 			\n" + 
	"uniform vec4  background; 		\n" + 
	"varying vec4  vColor; 			\n" + 
	"\n" 								+
	"void main() { 					\n" + 
	"   mat3  rotation_rhcp = mat3( 0.5, -0.5, 0.707, -0.5, 0.5, 0.707, -0.707, -0.707, 0.0); 	\n" + 
	"   //mat3  rotation_lhcp = mat3( 0.5, 0.5, 0.707, 0.5, 0.5, -0.707, -0.707, 0.707, 0.0); 		\n" + 
	"   mat3  rotation_lhcp = mat3( 0.5, 0.5, -0.707, 0.5, 0.5, 0.707, 0.707, -0.707, 0.0); 		\n" + 
	"	vec3 newPosition = position.xyz; \n" + 
	"	if(rotate) { 				\n" + 
	"		if(rhcp) { 				\n" + 
	"			newPosition = rotation_rhcp * position.xyz; \n" + 
	"		} else { 				\n" + 
	"			newPosition = rotation_lhcp * position.xyz; \n" + 
	"		} 						\n" + 
	"	} else { 					\n" + 
	"		newPosition = position.xyz; \n" + 
	"	} \n" + 
	"	float bright = mod(uTime + (swapAnimDir ? uv.y : uv.x), 1.0/wavesPerRing) * wavesPerRing; 			\n" + 
	"	if(animate) { 				\n" + 
	"		float shift = mod((uTime + uv.x), 1.0); 														\n" +
	"		if(shift <= cutoff) { \n" + 
	"			vColor = vec4( (1.0-bright)*color[0]+bright*background[0], (1.0-bright)*color[1]+bright*background[1], (1.0-bright)*color[2]+bright*background[2], color[3] );	\n" + 
	"		} else { \n" + 
	"			vColor = vec4( 0.2*color[0]+0.8*background[0], 0.2*color[1]+0.8*background[1], 0.2*color[2]+0.8*background[2], color[3] );	\n" + 
	"		}  \n" + 
	"		newPosition.z = newPosition.z - shift * 2.0 * radius;	\n" + 
	"	} else { 																							\n" + 
	"		float shift = mod((uTime + uv.x), 1.0); 														\n" +
	"		if(shift <= cutoff) { \n" + 
	"			//vColor = vec4( (1.0-bright)*color[0], (1.0-bright)*color[1], (1.0-bright)*color[2], color[3] );	\n" + 
	"			vColor = vec4( (1.0-bright)*color[0]+bright*background[0], (1.0-bright)*color[1]+bright*background[1], (1.0-bright)*color[2]+bright*background[2], color[3] );	\n" + 
	"		} else { \n" + 
	"			vColor = vec4( 0.0, 0.0, 0.0, 0.0 );	\n" + 
	"			newPosition = vec3( 0.0, 0.0, 0.0 );	\n" + 
	"		} 																								\n" + 
	"	} 																									\n" + 
	"	\n" 																									+	
	"	vec4 mvPosition = modelViewMatrix * vec4( newPosition.xyz*scale, 1.0 ); 							\n" + 
	"	gl_PointSize = 5.0; 																				\n" + 
	"	gl_Position = projectionMatrix * mvPosition; 														\n" + 
	"} \n";

// Heatmap vertex shader - Use this to update the position due to camera projection 
var densitymapVertexShader = 
	"uniform vec4  hiColor; 					\n" + 
	"uniform vec4  loColor;						\n" + 
	"uniform float maxHistogramValue;			\n" + 
	"uniform float contrast;					\n" + 
	"uniform float shape;						\n" + 
	"varying vec4 vColor;						\n" + 
	"\n" 											+ 
	"void main() {								\n" + 
	"	vec4 mvPosition;						\n" + 
	"	// UV.x contains histogram value		\n" + 
	"	// UV.y is identical to UV.x 			\n" + 
	"	if(maxHistogramValue > 0.0) {			\n" + 
	"		vColor = pow((uv.x / maxHistogramValue), contrast) * (hiColor-loColor) + loColor;	\n" + 
	"		mvPosition = modelViewMatrix * vec4( position.xyz*pow((uv.x / maxHistogramValue),shape), 1.0 );	\n" + 
	"	} else {								\n" + 
	"		vColor = hiColor; 					\n" + 
	"		mvPosition = modelViewMatrix * vec4( position.xyz, 1.0 ); 	\n" + 
	"	}										\n" + 
	"											\n" + 
	"	gl_PointSize = 1.0;						\n" + 
	"	gl_Position = projectionMatrix * mvPosition; \n" + 
	"}	\n";

// CVF fragment shader 
var cvfFragmentShader = 
	"varying vec4 vColor; 			\n" +
	"void main() { 					\n" +
	"	gl_FragColor = vColor; 		\n" +
	"} 	\n";

// Grid vertex shader - Use this to update the position due to camera projection 
var gridVertexShader = 
	"uniform float gridBrightness; 	\n" + 
	"uniform vec4 background; 		\n" + 
	"varying vec4 vColor; 			\n" + 

	"void main() { 					\n" + 
	"	//vColor = vec4( gridBrightness, gridBrightness, gridBrightness, 1.0 ); \n" + 
	"		vColor = vec4( (background[0]<0.5) ? gridBrightness*(1.0-background[0])+background[0] : background[0] - gridBrightness*background[0], " + 
	" 			(background[1]<0.5) ? gridBrightness*(1.0-background[1])+background[1] : background[1] - gridBrightness*background[1], " + 
	" 			(background[2]<0.5) ? gridBrightness*(1.0-background[2])+background[2] : background[2] - gridBrightness*background[2], 1.0 ); \n" + 
	"	vec4 mvPosition = modelViewMatrix * vec4( position.xyz, 1.0 ); 				\n" + 
	"	gl_PointSize = 1.0; 														\n" + 
	"	gl_Position = projectionMatrix * mvPosition; 								\n" + 
	"}";

var newGridVertexShader = 
	"uniform float gridBrightness; 	\n" + 
	"uniform float gridSize; 		\n" + 
	"uniform float gridDepth; 		\n" + 
	"uniform float waveSpeed; 		\n" + 
	"uniform float uTime; 			\n" + 
	"uniform vec4 background; 		\n" + 
	"varying vec4 vColor; 			\n" + 

	"void main() { 					\n" + 
	"	vec3 newPosition = position.xyz; 											\n" + 
	"	//newPosition.z = position.z * gridSize - mod(uTime*waveSpeed, gridSize); 	\n" + 
	"	newPosition.z = position.z * gridSize - mod(uTime*0.15, gridSize); 			\n" + 
	"	if(newPosition.z > -gridSize) {  											\n" + 
	"   	float br = abs(newPosition.z / gridSize) * gridBrightness; 				\n" +
	"		//vColor = vec4( br, br, br, 1.0 ); 								\n" + 
	"		vColor = vec4( (background[0]<0.5) ? br*(1.0-background[0])+background[0] : background[0] - br*background[0], " + 
	" 			(background[1]<0.5) ? br*(1.0-background[1])+background[1] : background[1] - br*background[1], " + 
	" 			(background[2]<0.5) ? br*(1.0-background[2])+background[2] : background[2] - br*background[2], 1.0 ); \n" + 
	"	} else if(newPosition.z < -(gridSize * gridDepth)) {  						\n" + 
	"   	float br = (1.0 - (abs(newPosition.z - (-gridSize * gridDepth)) / gridSize)) * gridBrightness; \n" +
	"		//vColor = vec4( br : br, br, br, alphaGrid ); 							\n" + 
	"		vColor = vec4( (background[0]<0.5) ? br*(1.0-background[0])+background[0] : background[0] - br*background[0], " + 
	" 			(background[1]<0.5) ? br*(1.0-background[1])+background[1] : background[1] - br*background[1], " + 
	" 			(background[2]<0.5) ? br*(1.0-background[2])+background[2] : background[2] - br*background[2], 1.0 ); \n" + 
	"	} else {																	\n" + 
	"		// vColor = vec4( gridBrightness, gridBrightness, gridBrightness, alphaGrid ); \n" + 
	"   	float br = gridBrightness; 				\n" +
	"		vColor = vec4( (background[0]<0.5) ? br*(1.0-background[0])+background[0] : background[0] - br*background[0], " + 
	" 			(background[1]<0.5) ? br*(1.0-background[1])+background[1] : background[1] - br*background[1], " + 
	" 			(background[2]<0.5) ? br*(1.0-background[2])+background[2] : background[2] - br*background[2], 1.0 ); \n" + 
	"	} 																			\n" + 
	"	newPosition.z += gridSize; \n" +
	"	vec4 mvPosition = modelViewMatrix * vec4( newPosition.xyz, 1.0 ); 			\n" + 
	"	gl_PointSize = 1.0; 														\n" + 
	"	gl_Position = projectionMatrix * mvPosition; 								\n" + 
	"}";


// Geometry object based on GUTCP equations, selected using "mode" variable:
//     BECVF (Eq 1.84)  : mode == 1 
//     OCVF  (Eq 1.95)  : mode == 2 
// Y00 BECVF (Eq 1.103) : mode == 3
// Y00 OCVF  (Eq 1.109) : mode == 4
class CVF {
	// Return an individual CVF vertex:
	// mode    = 1:BECVF, 2:OCVF
	// radius  = CVF radius (typically 100.0)
	// i_theta = iterator for theta
	// THETA   = number of current rings
	// i_phi   = iterator for PHI
	// PHI     = number of point in the current loop
	static getCvfVertex(mode, radius, i_theta, THETA, i_phi, PHI) {
		var theta = 2.0 * Math.PI * i_theta / THETA;

		var cvfRot = [];
		if(mode == 1) {
			// Setup the 3x3 rotation matrix from GUTCP Eq(1.84) BECVF
			cvfRot[0+0] = ( 0.5 + 0.5*Math.cos(theta));
			cvfRot[0+1] = (-0.5 + 0.5*Math.cos(theta));
			cvfRot[0+2] = (-0.70711*Math.sin(theta));
			cvfRot[3+0] = (-0.5 + 0.5*Math.cos(theta));
			cvfRot[3+1] = ( 0.5 + 0.5*Math.cos(theta)); 
			cvfRot[3+2] = (-0.70711*Math.sin(theta));
			cvfRot[6+0] = ( 0.70711*Math.sin(theta));
			cvfRot[6+1] = ( 0.70711*Math.sin(theta));
			cvfRot[6+2] = (Math.cos(theta));
		} else {
			// Setup the 3x3 rotation matrix from GUTCP Eq(1.95) OCVF
			cvfRot[0+0] = ( 0.25 * (1.0 + 3.0*Math.cos(theta)));
			cvfRot[0+1] = ( 0.25 * (-1.0 + Math.cos(theta) + 2.0 * 1.414 * Math.sin(theta)));
			cvfRot[0+2] = ( 0.25 * (-1.414 + 1.414 * Math.cos(theta) - 2.0 * Math.sin(theta)));
			cvfRot[3+0] = ( 0.25 * (-1.0 + Math.cos(theta) - 2.0 * 1.414 * Math.sin(theta)));
			cvfRot[3+1] = ( 0.25 * (1.0 + 3.0*Math.cos(theta))); 
			cvfRot[3+2] = ( 0.25 * (1.414 - 1.414 * Math.cos(theta) - 2.0 * Math.sin(theta)));
			cvfRot[6+0] = ( 0.5 * (((-1 + Math.cos(theta))/1.414) + Math.sin(theta)));
			cvfRot[6+1] = ( 0.25 * (1.414 - 1.414 * Math.cos(theta) + 2.0 * Math.sin(theta)));
			cvfRot[6+2] = (Math.cos(theta*0.5) * Math.cos(theta*0.5));
		}
		
		var phi;
		var cvf = [];
		if(mode == 1) {
			// Setup the BECVF basis current loop:
			phi = 2.0 * Math.PI * i_phi / PHI;
			cvf = [0.0, (radius * Math.cos(phi)), (-radius * Math.sin(phi))];
		} else {
			// Setup the OCVF basis current loop:
			phi = 2.0 * Math.PI * i_phi / PHI; // + 0.5 * Math.PI;  // Half-PI phase addition is to correct a kink visual artefact
			cvf = [(0.707*radius * Math.cos(phi)), 
				   (0.707*radius * Math.cos(phi)), 
				   (-radius * Math.sin(phi))];
		}
	
		var vertex = [];
		vertex[ 0 ] = cvf[0]*cvfRot[0+0] + cvf[1]*cvfRot[0+1] + cvf[2]*cvfRot[0+2];
		vertex[ 1 ] = cvf[0]*cvfRot[3+0] + cvf[1]*cvfRot[3+1] + cvf[2]*cvfRot[3+2]; 
		vertex[ 2 ] = cvf[0]*cvfRot[6+0] + cvf[1]*cvfRot[6+1] + cvf[2]*cvfRot[6+2];
	
		return vertex;
	}

	// Return an individual CVF vertex:
	// mode    = 1:BECVF, 2:OCVF
	// radius  = CVF radius (typically 100.0)
	// i_n     = iterator for N
	// N       = number of N rotations
	// i_m     = iterator for M
	// M       = number of M rotations
	// i_phi   = iterator for PHI
	// PHI     = number of point in the current loop
	static getY00Vertex(mode, radius, i_n, N, i_m, M, i_phi, PHI) {

		var m2piONm = 2.0 * Math.PI * i_m / M;
		var m_matrix = [];

		// The m_matrix is identical for both Y00 geometry types:
		m_matrix[0+0] = (0.25 * (1.0 + 3.0 * Math.cos(m2piONm)));
		m_matrix[0+1] = (0.25 * (-1.0 + Math.cos(m2piONm) + 2.828 * Math.sin(m2piONm)));
		m_matrix[0+2] = (0.25 * (-1.414 + 1.414 * Math.cos(m2piONm) - 2.0 * Math.sin(m2piONm)));
		m_matrix[3+0] = (0.25 * (-1.0 + Math.cos(m2piONm) - 2.828 * Math.sin(m2piONm)));
		m_matrix[3+1] = (0.25 * (1.0 + 3.0 * Math.cos(m2piONm)));
		m_matrix[3+2] = (0.25 * (1.414 - 1.414 * Math.cos(m2piONm) - 2.0 * Math.sin(m2piONm)));
		m_matrix[6+0] = (0.5 * (0.707 * (-1.0 + Math.cos(m2piONm)) + Math.sin(m2piONm)));
		m_matrix[6+1] = (0.25 * (1.414 - 1.414 * Math.cos(m2piONm) + 2.0 * Math.sin(m2piONm)));
		m_matrix[6+2] = (Math.cos(0.5 * m2piONm) * Math.cos(0.5 * m2piONm));
		
		var n2piONn = 2.0 * Math.PI * i_n / N;
		var n_matrix = [];
		// n_matrix is different, dependent on Y00 geometry type:
		if(mode == 1) {
				n_matrix[0+0] = (0.5 + 0.5 * Math.cos(n2piONn));
				n_matrix[0+1] = (-0.5 + 0.5 * Math.cos(n2piONn));
				n_matrix[0+2] = (-Math.sin(n2piONn) * 0.707);
				n_matrix[3+0] = (-0.5 + 0.5 * Math.cos(n2piONn));
				n_matrix[3+1] = (0.5 + 0.5 * Math.cos(n2piONn));
				n_matrix[3+2] = (-Math.sin(n2piONn) * 0.707);
				n_matrix[6+0] = (Math.sin(n2piONn) * 0.707);
				n_matrix[6+1] = (Math.sin(n2piONn) * 0.707);
				n_matrix[6+2] = (Math.cos(n2piONn));
		} else {
				n_matrix[0+0] = (0.5 * Math.cos(n2piONn) - 0.5 * Math.sin(n2piONn));
				n_matrix[0+1] = (0.707 * Math.cos(n2piONn) + 0.707 * Math.sin(n2piONn));
				n_matrix[0+2] = (0.5 * Math.cos(n2piONn) - 0.5 * Math.sin(n2piONn));
				n_matrix[3+0] = (-0.5 * Math.cos(n2piONn) - 0.5 * Math.sin(n2piONn));
				n_matrix[3+1] = (0.707 * Math.cos(n2piONn) - 0.707 * Math.sin(n2piONn));
				n_matrix[3+2] = (-0.5 * Math.cos(n2piONn) - 0.5 * Math.sin(n2piONn));
				n_matrix[6+0] = (-0.707);
				n_matrix[6+1] = (0.0);
				n_matrix[6+2] = (0.707);
		}

		var phi = 2.0 * Math.PI * i_phi / PHI;
		var becvf ;
		if(mode == 1) {
			becvf = [(0.0),
					 (radius * Math.cos(phi)),
					 (-radius * Math.sin(phi))];
		} else {
			becvf = [(radius * Math.cos(phi)),
					 (radius * Math.sin(phi)),
					 (0.0)];
		}
		var tx = n_matrix[0+0] * becvf[0] + n_matrix[0+1] * becvf[1] + n_matrix[0+2] * becvf[2];
		var ty = n_matrix[3+0] * becvf[0] + n_matrix[3+1] * becvf[1] + n_matrix[3+2] * becvf[2];
		var tz = n_matrix[6+0] * becvf[0] + n_matrix[6+1] * becvf[1] + n_matrix[6+2] * becvf[2];

		var vertex = [];
		vertex[0] = m_matrix[0+0] * tx + m_matrix[0+1] * ty + m_matrix[0+2] * tz;
		vertex[1] = m_matrix[3+0] * tx + m_matrix[3+1] * ty + m_matrix[3+2] * tz;
		vertex[2] = m_matrix[6+0] * tx + m_matrix[6+1] * ty + m_matrix[6+2] * tz;
		return vertex;
	}
	
	// Return an individual CVF vertex:
	// mode    = 1:RH, 2:LH
	// field   = 1:Electric, 2:Magnetic
	// radius  = CVF radius (typically 100.0)
	// i_theta = iterator for theta
	// THETA   = number of current rings
	// i_phi   = iterator for PHI
	// PHI     = number of point in the current loop
	static getPhotonVertex(mode, field, radius, i_theta, THETA, i_phi, PHI) {
		//var theta = 2.0 * Math.PI * i_theta / THETA;
		var cos_theta = Math.cos(2.0 * Math.PI * i_theta / THETA);
		var sin_theta = Math.sin(2.0 * Math.PI * i_theta / THETA);

		var cvfRot = [];
		if(mode == 1) {
			// Setup the 3x3 rotation matrix from GUTCP Eq(V-4) Right-Handed
			cvfRot[0+0] = ( 0.5 + 0.5*cos_theta);
			cvfRot[0+1] = ( 0.5 - 0.5*cos_theta);
			cvfRot[0+2] = (-0.70711*sin_theta);
			cvfRot[3+0] = ( 0.5 - 0.5*cos_theta);
			cvfRot[3+1] = ( 0.5 + 0.5*cos_theta); 
			cvfRot[3+2] = ( 0.70711*sin_theta);
			cvfRot[6+0] = ( 0.70711*sin_theta);
			cvfRot[6+1] = (-0.70711*sin_theta);
			cvfRot[6+2] = (cos_theta);
		} else {
			// Setup the 3x3 rotation matrix from GUTCP Eq(V-9) Left-Handed
			cvfRot[0+0] = ( 0.5 + 0.5*cos_theta);
			cvfRot[0+1] = (-0.5 + 0.5*cos_theta);
			cvfRot[0+2] = ( 0.70711*sin_theta);
			cvfRot[3+0] = (-0.5 + 0.5*cos_theta);
			cvfRot[3+1] = ( 0.5 + 0.5*cos_theta); 
			cvfRot[3+2] = ( 0.70711*sin_theta);
			cvfRot[6+0] = (-0.70711*sin_theta);
			cvfRot[6+1] = (-0.70711*sin_theta);
			cvfRot[6+2] = (cos_theta);
		}
		
		var phi;
		var cvf = [];
		if(field == 1) {
			// Setup the electric field loop:
			phi = 2.0 * Math.PI * i_phi / PHI;
			cvf = [0.0, (radius * Math.cos(phi)), (radius * Math.sin(phi))];
		} else {
			// Setup the magnetic field loop:
			phi = 2.0 * Math.PI * i_phi / PHI; // + 0.5 * Math.PI;  // Half-PI phase addition is to correct a kink visual artefact
			cvf = [(radius * Math.cos(phi)), 0.0, (radius * Math.sin(phi))];
		}
	
		var vertex = [];
		vertex[ 0 ] = cvf[0]*cvfRot[0+0] + cvf[1]*cvfRot[0+1] + cvf[2]*cvfRot[0+2];
		vertex[ 1 ] = cvf[0]*cvfRot[3+0] + cvf[1]*cvfRot[3+1] + cvf[2]*cvfRot[3+2]; 
		vertex[ 2 ] = cvf[0]*cvfRot[6+0] + cvf[1]*cvfRot[6+1] + cvf[2]*cvfRot[6+2];
		
		return vertex;
		
		// Rotation of the photon fields:
		var vertex2 = [];
		if(mode == 1) {
			// Now rotate RH about (-1, 1, 0) axis by -90 degrees, so that it can move as a TEM wave in the +z direction:
			// [Or equivalently, about (1, -1, 0) by 90 degrees.]
			vertex2[ 0 ] =    0.5*vertex[0] -   0.5*vertex[1] - 0.707*vertex[2];
			vertex2[ 1 ] =   -0.5*vertex[0] +   0.5*vertex[1] - 0.707*vertex[2]; 
			vertex2[ 2 ] =  0.707*vertex[0] + 0.707*vertex[1] + 0.0;
		} else {
			// Now rotate LH about (1, 1, 0) axis by +90 degrees, so that it can move as a TEM wave in the +z direction:
			vertex2[ 0 ] =    0.5*vertex[0] +   0.5*vertex[1] + 0.707*vertex[2];
			vertex2[ 1 ] =    0.5*vertex[0] +   0.5*vertex[1] - 0.707*vertex[2]; 
			vertex2[ 2 ] = -0.707*vertex[0] + 0.707*vertex[1] + 0.0;
		}
		return vertex2;
	}

	// Returns a two-element array containing geometry and uv for GUTCP Current-Vector Fields.
	// mode 1 - returns a BECVF geometries and uvs
	// mode 2 - returns a OCVF geometry and uvs 
	static createCvfGeometry (mode, radius, THETA, PHI) {
		// Storage for all the locations for the geometry:
		var positions;
		var uvs;
	
		// Local loop variables:
		var retval = [];
		retval[0] = [];
		retval[1] = [];
		var vtx = [];
		for ( var i = 0; i < THETA; i++ ) {
			var p = 0;
			var q = 0;
			positions = new Float32Array( PHI * 3 );
			uvs = new Float32Array( PHI * 2 );
			for ( var j = 0; j < PHI; j++ ) {
				vtx = CVF.getCvfVertex(mode, radius, i, THETA, j, PHI);
				positions[ p++ ] = vtx[0];
				positions[ p++ ] = vtx[1];
				positions[ p++ ] = vtx[2];
				
				// Setup the UV array, which will contain a lookup-table of 0.0-1.0 for (i,j) based on the index:
				uvs[ q++ ] = i / ( THETA - 1 );
				uvs[ q++ ] = j / ( PHI - 1 );
			}
			retval[0].push(positions);
			retval[1].push(uvs);
		}
		return retval;
	}

	// Static function that returns a two-element array containing geometry and uv for GUTCP Current-Vector Fields.
	// Y00 BECVF (Eq 1.103) : mode == 1
	// Y00 OCVF  (Eq 1.109) : mode == 2
	static createY00Geometry (mode, radius, N, M, PHI) {
		// Storage for all the locations for the BECVF:
		var positions;
		var uvs;
		
		var retval = [];
		retval[0] = [];
		retval[1] = [];
		var vtx = [];
		for (var i_m = 1; i_m <= M; i_m++) {
			for (var i_n = 1; i_n <= N; i_n++) {
				var p = 0;
				var q = 0;
				positions = new Float32Array( PHI * 3 );
				uvs = new Float32Array( PHI * 2 );
				for (var i_phi = 0; i_phi < PHI; i_phi++) {
					vtx = CVF.getY00Vertex(mode, radius, i_n, N, i_m, M, i_phi, PHI);
					positions[p++] = vtx[0];
					positions[p++] = vtx[1];
					positions[p++] = vtx[2];
					
					uvs[ q++ ] = i_phi / ( PHI - 1 );
					uvs[ q++ ] = i_phi / ( PHI - 1 );
				}
				retval[0].push(positions);
				retval[1].push(uvs);
			}
		}
		return retval;
	}
	
	// Static function that returns a two-element array containing geometry and uv for GUTCP Current-Vector Fields.
	// Y00 BECVF (Eq 1.103) : mode == 1
	// Y00 OCVF  (Eq 1.109) : mode == 2
	static createY00Geometry (mode, radius, N, M, PHI) {
		// Storage for all the locations for the BECVF:
		var positions;
		var uvs;
		
		var retval = [];
		retval[0] = [];
		retval[1] = [];
		var vtx = [];
		for (var i_m = 1; i_m <= M; i_m++) {
			for (var i_n = 1; i_n <= N; i_n++) {
				var p = 0;
				var q = 0;
				positions = new Float32Array( PHI * 3 );
				uvs = new Float32Array( PHI * 2 );
				for (var i_phi = 0; i_phi < PHI; i_phi++) {
					vtx = CVF.getY00Vertex(mode, radius, i_n, N, i_m, M, i_phi, PHI);
					positions[p++] = vtx[0];
					positions[p++] = vtx[1];
					positions[p++] = vtx[2];
					
					uvs[ q++ ] = i_phi / ( PHI - 1 );
					uvs[ q++ ] = i_phi / ( PHI - 1 );
				}
				retval[0].push(positions);
				retval[1].push(uvs);
			}
		}
		return retval;
	}
	
	// Static function that returns a two-element array containing geometry and uv for GUTCP Current-Vector Fields.
	// Photon Right Handed : mode == 1
	// Photon Left Handed  (Eq 1.109) : mode == 2
	static createPhotonGeometry (mode, field, radius, THETA, PHI) {
		// Storage for all the locations for the geometry:
		var positions;
		var uvs;
	
		// Local loop variables:
		var retval = [];
		retval[0] = [];
		retval[1] = [];
		var vtx = [];
		for ( var i = 0; i < THETA; i++ ) {
			var p = 0;
			var q = 0;
			positions = new Float32Array( PHI * 3 );
			uvs = new Float32Array( PHI * 2 );
			for ( var j = 0; j < PHI; j++ ) {
				vtx = CVF.getPhotonVertex(mode, field, radius, i, THETA, j, PHI);
				positions[ p++ ] = vtx[0];
				positions[ p++ ] = vtx[1];
				positions[ p++ ] = vtx[2];
				
				// Setup the UV array, which will contain a lookup-table of 0.0-1.0 for (i,j) based on the index:
				uvs[ q++ ] = i / ( THETA - 1 );
				uvs[ q++ ] = j / ( PHI - 1 );
			}
			retval[0].push(positions);
			retval[1].push(uvs);
		}
		return retval;
	}
	
	// Construct a BECVF object with a radius of "radius", a rotation resolution of "THETA", and a 
	// great-circle current-loop density of "PHI":
	constructor(mode, radius, THETA, PHI) {
		this.scene = 0;
		// Variable handles for dat.gui to manipulate:
		this.color = [255, 0, 0, 1.0]; 	// Let's default to RED in [R, G, B, A]
		this.animate = true;
		this.swapAnimDir = false;
		this.visibility = true;
		this.scale = 1.0;
		this.sphere = 1.0;
		this.pole = [0.0, -radius, 0.0];
		this.cutoff = 0.25;
		this.radius = radius;
		this.background = [0.0, 0.0, 0.0];
		this.rotate = false;
		this.rhcp = false;
		
		// Create the shader uniforms:
		this.cvfUniforms = {
			//cameraConstant: { value: getCameraConstant( camera ) },
			uTime:     	{ value: 0.0 },
			waveSpeed: 	{ value: 0.25 },
			wavesPerRing: { value: 1.0 },
			background: { value: new THREE.Vector4(this.background[0]/255.0, this.background[1]/255.0, this.background[2]/255.0, 1.0) },
			radius:  	{ value: this.radius },
			scale:   	{ value: this.scale },
			sphere:  	{ value: this.sphere },
			pole:  	 	{ value: this.pole },
			cutoff:  	{ value: this.cutoff },
			animate: 	{ value: this.animate },
			swapAnimDir: { value: this.swapAnimDir},
			rotate: 	{ value: this.rotate },
			rhcp: 		{ value: this.rhcp },
			color:   	{ value: new THREE.Vector4(this.color[0]/255.0, this.color[1]/255.0, this.color[2]/255.0, this.color[3]) }
		};
		
		this.theta = THETA;
		this.N_M = THETA;
		this.phi = PHI;
		this.loadNewGeometry(mode, radius, THETA, PHI);
	}
	
	// Insert the geometries into the scenegraph's sceneObject:
	insertScene(sceneObject) {
		sceneObject.add(this.cvfLines);
		this.scene = sceneObject;
		if(this.scene.background) {
			this.background = this.scene.background;
		} 
	}
	
	// Loads a new geometry based on the given set of parameters. Hopefully previous geometries are unloaded from
	// the GPU's VRAM.
	loadNewGeometry(mode, radius, THETA, PHI) {
		// Remove the geometries from the scene:
		if(this.scene) {
			this.scene.remove(this.cvfLines);
		}
		
		var vertex_shader = cvfVertexShader;
		
		// Re-create the geometries:
		var cvfGeo;
		if((mode == 1) || (mode == 2)) {
			cvfGeo = CVF.createCvfGeometry (mode, radius, THETA, PHI);
		} else 
		// assume mode == 3 or 4 if we get here
		if((mode == 3) || (mode == 4)) {
			cvfGeo = CVF.createY00Geometry (mode-2, radius, THETA/5, THETA/5, PHI);
		} else 
		if(mode == 5) {
			cvfGeo = CVF.createPhotonGeometry (1, 1, radius, THETA, PHI);
			vertex_shader = photonVertexShader;
			this.cvfUniforms.rhcp.value = true;
		} else 
		if(mode == 6) {
			cvfGeo = CVF.createPhotonGeometry (1, 2, radius, THETA, PHI);
			vertex_shader = photonVertexShader;
			this.cvfUniforms.rhcp.value = true;
		} else 
		if(mode == 7) {
			cvfGeo = CVF.createPhotonGeometry (2, 1, radius, THETA, PHI);
			vertex_shader = photonVertexShader;
			this.cvfUniforms.rhcp.value = false;
		} else 
		if(mode == 8) {
			cvfGeo = CVF.createPhotonGeometry (2, 2, radius, THETA, PHI);
			vertex_shader = photonVertexShader;
			this.cvfUniforms.rhcp.value = false;
		} else { // ERROR condition!!!
			cvfGeo = CVF.createCvfGeometry (1, radius, THETA, PHI);
		} 
		
		// The BECVF geometry, based on BufferGeometry for efficiency:
		this.positions = cvfGeo[0];
		this.uvs = cvfGeo[1];

		// Create and initialise the ShaderMaterial:
		this.material = new THREE.ShaderMaterial( {
			uniforms:       this.cvfUniforms,
			vertexShader:   vertex_shader,
			fragmentShader: cvfFragmentShader
		} );
		this.material.extensions.drawBuffers = true;

		this.cvfLines = new THREE.Object3D();
		for(var p=0, pl=this.positions.length; p<pl; ++p) {
			this.geometry = new THREE.BufferGeometry();
			// Add the vertices and uv arrays as attributes into the geometry:
			this.geometry.addAttribute( 'position', new THREE.BufferAttribute( this.positions[p], 3 ) );
			this.geometry.addAttribute( 'uv', new THREE.BufferAttribute( this.uvs[p], 2 ) );
			
			// Create the Lines geometry and load the into the scene:
			var lines = new THREE.LineLoop( this.geometry, this.material );
			lines.matrixAutoUpdate = false;
			lines.updateMatrix();
			this.cvfLines.add(lines);
		}
		
		this.setVisibility(this.visibility);
		
		// Add the new geometries:
		if(this.scene) {
			this.scene.add(this.cvfLines);
		}
	}
	
	// Set the visible geometry. No point ever displaying both, as they use the SAME vertices:
	setVisibility(value) {
		this.visibility = value;
		this.cvfLines.visible = value;
	}
	
	// Select whether to animate or not:
	setAnimate(value) {
		this.cvfUniforms.animate.value = (value) ? true : false;
	}
	
	// Select whether to animate or not:
	setRotate(value) {
		this.cvfUniforms.rotate.value = (value) ? true : false;
	}
	
	setSpeed(value) {
		this.cvfUniforms.waveSpeed.value = value;
	}
	
	// Little hack for dat.gui. Update the color by passing 3-element array in form [255, 128, 0]:
	setColor(value) {
		this.color = value;
		this.cvfUniforms.color.value = [this.color[0]/255.0, this.color[1]/255.0, this.color[2]/255.0, 1.0];
	}
	
	// Little hack for dat.gui. Update the color by passing 3-element array in form [255, 128, 0]:
	setBackground(value) {
		this.background = value;
		this.cvfUniforms.background.value = [this.background[0]/255.0, this.background[1]/255.0, this.background[2]/255.0, 1.0];
	}
	
	// Update the time for the shader to see:
	addTime( delta_time ) {
		this.cvfUniforms.uTime.value += delta_time * this.cvfUniforms.waveSpeed.value;
	}
	
	// Update the time for the shader to see:
	setTime( tick_time ) {
		this.cvfUniforms.uTime.value = tick_time;
	}
	
	swapAnimationDirection(value) {
		this.cvfUniforms.swapAnimDir.value = value;
	}
}

// Density map class for calculating the CVF heat/density map:
class DensityMap {
	constructor() {
		// Heatmap specific stuff:
		this.densitymapGeometry  = new THREE.IcosahedronBufferGeometry( 98, 5 );
		this.densitymapVertices  = this.densitymapGeometry.getAttribute('position');
		this.densitymapHistogram = this.densitymapGeometry.getAttribute('uv');
		//var densitymapColor = new Float32Array( densitymapVertices.count * 3 );
	
		// Initialise to ZERO the UV values:
		for(var i = 0; i < this.densitymapHistogram.count; i++) {
			this.densitymapHistogram.setX(i, 0.0);
			this.densitymapHistogram.setY(i, 0.0);
		}

		this.hiColor = [0, 255, 0, 1.0]; 	// Let's default to RED in [R, G, B, A]
		this.loColor = [0, 0, 255, 1.0]; 	// Let's default to RED in [R, G, B, A]
		this.contrast = 1.0;
		this.shape = 1.0;
		this.densitymapUniforms = {
			cameraConstant: 	{ value: getCameraConstant( camera ) },
			hiColor: 			{ value: new THREE.Vector4(this.hiColor[0]/255.0, this.hiColor[1]/255.0, this.hiColor[2]/255.0, this.hiColor[3]) },
			loColor: 			{ value: new THREE.Vector4(this.loColor[0]/255.0, this.loColor[1]/255.0, this.loColor[2]/255.0, this.loColor[3]) },
			maxHistogramValue: 	{ value: 0.0 },
			contrast:			{ value: this.contrast },
			shape: 				{ value: this.shape }
		};

		// ShaderMaterial
		this.wireframe = true;
		this.material = new THREE.ShaderMaterial( {
												uniforms:       this.densitymapUniforms,
												vertexShader:   densitymapVertexShader,
												fragmentShader: cvfFragmentShader,
												wireframe: 		this.wireframe
												} );
		this.material.extensions.drawBuffers = true;
		this.densitymapMesh = new THREE.Mesh( this.densitymapGeometry, this.material );
		this.visibility = true;
		this.densitymapMesh.visible = this.visibility;
		
		this.N = 12;
		this.M = 12;
		this.PHI = 60;
		this.hasWork = false;
		
		// The next two for-loops are to create a 3-dimensional lookup-table to significantly
		// speed-up the generation of the heat-map by only iterating through a limited set of
		// vertices. Without this, high-resolution spheres become unwieldy, as every vertex in the 
		// sphere has to be ranged for EVERY point in the NxMxPHI orbitsphere. (Poor scalability.)
		
		// Initialise the 12x12x12 lookup table indices:
		this.fast_index = [];
		for(var i_x=0; i_x < 12; i_x++) {
			this.fast_index[i_x] = [];
			for(var i_y=0; i_y < 12; i_y++) {
				this.fast_index[i_x][i_y] = [];
				for(var i_z=0; i_z < 12; i_z++) {
					this.fast_index[i_x][i_y][i_z] = [];
				}
			}
		}
		// Populate the lookup table with the vertex and its index:
		for(var i = 0; i < this.densitymapVertices.count; i++) {
			var hx = Math.trunc(this.densitymapVertices.getX(i)/25.0) + 6;
			var hy = Math.trunc(this.densitymapVertices.getY(i)/25.0) + 6;
			var hz = Math.trunc(this.densitymapVertices.getZ(i)/25.0) + 6;
			this.fast_index[hx][hy][hz].push(i);
		}
	}
	
	// Insert the geometries into the scenegraph's sceneObject:
	insertScene(sceneObject) {
		sceneObject.add(this.densitymapMesh);
	}
	
	// Set the visible geometry. No point ever displaying both, as they use the SAME vertices:
	setVisibility(value) {
		this.visibility = value;
		this.densitymapMesh.visible = value;
	}
	
	// Select whether to animate or not:
	setAnimate(value) {
		this.densitymapUniforms.animate.value = (value) ? true : false;
	}
	
	// Little hack for dat.gui. Update the color by passing 4-element array in form [255, 128, 0, 1.0]:
	setHiColor(value) {
		this.hiColor = value;
		this.densitymapUniforms.hiColor.value = [this.hiColor[0]/255.0, this.hiColor[1]/255.0, this.hiColor[2]/255.0, this.hiColor[3]];
	}
	
	setLoColor(value) {
		this.loColor = value;
		this.densitymapUniforms.loColor.value = [this.loColor[0]/255.0, this.loColor[1]/255.0, this.loColor[2]/255.0, this.loColor[3]];
	}
	
	reset() {
		// Clear the geometry UV:
		for(var i = 0; i < this.densitymapHistogram.count; i++) {
			this.densitymapHistogram.setX(i, 0.0);
			this.densitymapHistogram.setY(i, 0.0);
		}
		this.densitymapUniforms.maxHistogramValue.value = 0.0;
		this.densitymapHistogram.needsUpdate = true;
		this.i_n = 0;
		this.i_m = 0;
		this.i_phi = 0;
		this.hasWork = true;
	}
	
	pause() {
		this.hasWork = !this.hasWork;
	}
	
	// Function that solves the density map, but is capable of solving it in chunks of work, so as not to stall
	// the screen refresh cycle. That way, high values of N x M x Phi, which can take minutes to solve, but 
	// still allows a dynamic and responsive user interface.
	process(value) {
		if(this.hasWork) {
			var count = 0;
			OUT:
			for(; this.i_phi < this.PHI; ++this.i_phi) {
				for(; this.i_m < this.M; ++this.i_m) {
					for(; this.i_n < this.N; ++this.i_n) {
						var geo = CVF.getY00Vertex (1, 98.0, this.i_n, this.N, this.i_m, this.M, this.i_phi, this.PHI);

						// Obtain xyz coords of the CVF vertex:
						var hx = Math.trunc(geo[0]/25.0) + 6;
						var hy = Math.trunc(geo[1]/25.0) + 6;
						var hz = Math.trunc(geo[2]/25.0) + 6;
						
						// Look in the buckets/bins on either side of this point, to capture all vertices in range.
						// This means looking in a 3x3 grid of bins to search. Still much quicker then iterating 
						// through ALL vertices:
						for(var xx=-1; xx<2; ++xx) {
							for(var yy=-1; yy<2; ++yy) {
								for(var zz=-1; zz<2; ++zz) {
								
									// Now find all the vertices in the bin, and increment its histogram:
									for(var j=0, jl=this.fast_index[hx+xx][hy+yy][hz+zz].length; j < jl; ++j) {
										var pos = this.fast_index[hx+xx][hy+yy][hz+zz][j];
										var rangeSquared = Math.pow(geo[0] - this.densitymapVertices.getX(pos), 2)
														 + Math.pow(geo[1] - this.densitymapVertices.getY(pos), 2)
														 + Math.pow(geo[2] - this.densitymapVertices.getZ(pos), 2);
							
										// Use a range of 15 units, so 15*15=225:
										if(rangeSquared < 225.0) {
											var num = this.densitymapHistogram.getX(pos);
											this.densitymapHistogram.setX(pos, num+1.0);
		
											// Find highest value in array:
											if((num+1.0) > this.densitymapUniforms.maxHistogramValue.value) {
												this.densitymapUniforms.maxHistogramValue.value = num+1.0;
											}
										}
									}
									
								}
							}
						}
						
						// This is where we check if we've completed a certain number of searches. If we have
						// done more than "value", then we finish for now, and continue during the next
						// process() function call, and continue where we left-off:
						if(++count >= value) {
							this.densitymapHistogram.needsUpdate = true;
							break OUT;
						}
					}
					// Finished the N-loop, so re-zero the iterator here:
					this.i_n = 0;
				} 
				// Finished the M-loop, so re-zero the iterator here:
				this.i_m = 0;
			}
			// This checks the condition that we are completely done iterating through NxMxPHI. The needsUpdate
			// variable is how we tell THREE.js to update/reload the geometry in VRAM (GPU buffers)
			if(this.i_phi >= this.PHI) {
				this.hasWork = false;
				this.densitymapHistogram.needsUpdate = true;
			}
		}
	}
}

class Grid {
	constructor() {
		this.width = 200;
		this.height = 200;
		this.depth = 200;
		this.separators = 3;
		//this.generateBoundary();
		this.generateFull();
	}
	
	generateFull() {
		var spacing = this.width / this.separators; // pixels
		this.geometry = new THREE.BufferGeometry();

		// Storage for all the locations for the BECVF:
		var positions = new Float32Array( 6 * ((this.separators+1)**2) * 3 );
		var p = 0;
		var x = 0.0;
		var y = 0.0;
		var z = 0.0;

		// Populate the grid positions:
		for ( var i = 0; i <= this.separators; i++ ) {
			x = (i * spacing) - (this.separators * spacing * 0.5);
			for ( var j = 0; j <= this.separators; j++ ) {
				y = (j * spacing) - (this.separators * spacing * 0.5);

				positions[ p++ ] = x;
				positions[ p++ ] = y;
				positions[ p++ ] = -(this.separators * spacing * 0.5);

				positions[ p++ ] = x;
				positions[ p++ ] = y;
				positions[ p++ ] = (this.separators * spacing * 0.5);
			}
		}

		for ( var i = 0; i <= this.separators; i++ ) {
			x = (i * spacing) - (this.separators * spacing * 0.5);
			for ( var j = 0; j <= this.separators; j++ ) {
				y = (j * spacing) - (this.separators * spacing * 0.5);

				positions[ p++ ] = y;
				positions[ p++ ] = -(this.separators * spacing * 0.5);
				positions[ p++ ] = x;

				positions[ p++ ] = y;
				positions[ p++ ] = (this.separators * spacing * 0.5);
				positions[ p++ ] = x;
			}
		}

		for ( var i = 0; i <= this.separators; i++ ) {
			x = (i * spacing) - (this.separators * spacing * 0.5);
			for ( var j = 0; j <= this.separators; j++ ) {
				y = (j * spacing) - (this.separators * spacing * 0.5);

				positions[ p++ ] = -(this.separators * spacing * 0.5);
				positions[ p++ ] = x;
				positions[ p++ ] = y;

				positions[ p++ ] = (this.separators * spacing * 0.5);
				positions[ p++ ] = x;
				positions[ p++ ] = y;
			}
		}

		this.geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
		
		this.gridBrightness = 0.9;
		this.background = [0, 0, 0];
		
		this.gridUniforms = { 
			gridBrightness: { value: this.gridBrightness },
			background: 	{ value: new THREE.Vector4(this.background[0]/255.0, this.background[1]/255.0, this.background[2]/255.0, 1.0) }
		};
	
		// ShaderMaterial
		this.material = new THREE.ShaderMaterial( {
												uniforms:       this.gridUniforms,
												linewidth:      2,
												vertexShader:   gridVertexShader,
												fragmentShader: cvfFragmentShader
												} );                                                        
		this.material.extensions.drawBuffers = true;
	
		this.gridGeometry = new THREE.LineSegments( this.geometry, this.material );
		this.gridGeometry.matrixAutoUpdate = false;
		this.gridGeometry.updateMatrix();
	}
	
	generateBoundary() {
		// Create a two line-squares, followed by four lines after.
		var positions = new Float32Array( 24*3 ); // 12 boundary lines require 2 vertices each, 3 floats per vertex
		var x = 0.0;
		var y = 0.0;
		var p = 0;
		
		for ( var i = 0; i < 2; i++ ) {
			x = -this.width * 0.5 + i * this.width;
			for ( var j = 0; j < 2; j++ ) {
				y = -this.height * 0.5 + j * this.height;
				
				positions[p++] = x;
				positions[p++] = y;
				positions[p++] = -this.depth * 0.5;
		
				positions[p++] = x;
				positions[p++] = y;
				positions[p++] = this.depth * 0.5;
			}
		}
		
		for ( var i = 0; i < 2; i++ ) {
			x = -this.depth * 0.5 + i * this.depth;
			for ( var j = 0; j < 2; j++ ) {
				y = -this.height * 0.5 + j * this.height;
				
				positions[p++] = -this.width * 0.5;
				positions[p++] = y;
				positions[p++] = x;
		
				positions[p++] = this.width * 0.5;
				positions[p++] = y;
				positions[p++] = x;
			}
		}
		
		for ( var i = 0; i < 2; i++ ) {
			x = -this.width * 0.5 + i * this.width;
			for ( var j = 0; j < 2; j++ ) {
				y = -this.depth * 0.5 + j * this.depth;
				
				positions[p++] = x;
				positions[p++] = -this.height * 0.5;
				positions[p++] = y;
		
				positions[p++] = x;
				positions[p++] = this.height * 0.5;
				positions[p++] = y;
			}
		}
		
		this.geometry = new THREE.BufferGeometry();
		this.geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
		this.gridBrightness = 0.9;
		this.background = [0, 0, 0];
		
		this.gridUniforms = { 
			gridBrightness: { value: this.gridBrightness },
			background: 	{ value: new THREE.Vector4(this.background[0]/255.0, this.background[1]/255.0, this.background[2]/255.0, 1.0) }
		};
		
		// ShaderMaterial
		this.material = new THREE.ShaderMaterial( {
												uniforms:       this.gridUniforms,
												linewidth:      2,
												vertexShader:   gridVertexShader,
												fragmentShader: cvfFragmentShader
												} );                                                        
		this.material.extensions.drawBuffers = true;
		
		this.gridGeometry = new THREE.LineSegments( this.geometry, this.material );
		this.gridGeometry.matrixAutoUpdate = false;
		this.gridGeometry.updateMatrix();
	}
	
	setBackground(background) {
		this.gridUniforms.background.value = [background[0]/255.0, background[1]/255.0, background[2]/255.0, 1.0];
	}
	
	// Set the brightness of the grid. Make it invisible if brightness is 0.05 or less.
	setBrightness(brightness) {
		this.gridBrightness = brightness;
		this.gridUniforms.gridBrightness.value = brightness;
		if(brightness >= 0.05) {
			this.gridGeometry.visible = true;
		} else {
			this.gridGeometry.visible = false;
		}
	}
	
	// Insert the geometries into the scenegraph's sceneObject:
	insertScene(sceneObject) {
		sceneObject.add(this.gridGeometry);
	}
}

class NewGrid {
	constructor(width, height, depth, cube_size) {
		this.width = width; 		// in cubes
		this.height = height; 		// in cubes
		this.depth = depth; 		// in cubes
		this.cube_size = cube_size; // in pixels
		
		this.geometry = new THREE.BufferGeometry();

		// Storage for all the locations for the BECVF:
		var positions = new Float32Array( 6 * 2 * (this.width+1) * (this.height+1) * (this.depth+1));
		var p = 0;
		
		//var movers = [];
		
		// Vertical grid elements first:
		for(var z=0; z <= this.depth; ++z) {
			for(var y = 0; y <= this.height; ++y) {
				for(var x = 0; x <= this.width; ++x) {
					positions[ p++ ] = this.cube_size * (x - 0.5 * this.width);
					positions[ p++ ] = 0.5 * this.height * this.cube_size;
					positions[ p++ ] = -z;

					positions[ p++ ] = this.cube_size * (x - 0.5 * this.width);
					positions[ p++ ] = -0.5 * this.height * this.cube_size;
					positions[ p++ ] = -z;
				}
			}
		}
		
		// Now, do the horizontal grid elements next:
		for(var z=0; z <= this.depth; ++z) {
			for(var y = 0; y <= this.height; ++y) {
				for(var x = 0; x <= this.width; ++x) {
					positions[ p++ ] = -0.5 * this.width * this.cube_size;
					positions[ p++ ] = this.cube_size * (y - 0.5 * this.height);
					positions[ p++ ] = -z;

					positions[ p++ ] = 0.5 * this.width * this.cube_size;
					positions[ p++ ] = this.cube_size * (y - 0.5 * this.height);
					positions[ p++ ] = -z;
				}
			}
		}
		this.geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
		
		this.gridBrightness = 0.5;
		this.background = [0, 0, 0];
		
		this.gridUniforms = { 
			gridBrightness: { value: this.gridBrightness },
			gridSize: 		{ value: this.cube_size }, 
			gridDepth: 		{ value: this.depth }, 
			waveSpeed: 		{ value: 1.0 }, 
			alphaGrid: 		{ value: 1.0 },
			uTime: 			{ value: 0.0 },
			background: 	{ value: new THREE.Vector4(this.background[0]/255.0, this.background[1]/255.0, this.background[2]/255.0, 1.0) }
		};
	
		// ShaderMaterial
		this.material = new THREE.ShaderMaterial( {
												uniforms:       this.gridUniforms,
												linewidth:      2,
												vertexShader:   newGridVertexShader,
												fragmentShader: cvfFragmentShader
												} );                                                        
		this.material.extensions.drawBuffers = true;
	
		this.gridGeometry = new THREE.LineSegments( this.geometry, this.material );
		this.gridGeometry.matrixAutoUpdate = false;
		this.gridGeometry.updateMatrix();
	}
	
	insertScene(sceneObject) {
		sceneObject.add(this.gridGeometry);
	}
	
	setBackground(background) {
		this.gridUniforms.background.value = [background[0]/255.0, background[1]/255.0, background[2]/255.0, 1.0];
	}
	
	setBrightness(brightness) {
		this.gridBrightness = brightness;
		this.gridUniforms.gridBrightness.value = brightness;
		if(brightness >= 0.05) {
			this.gridGeometry.visible = true;
		} else {
			this.gridGeometry.visible = false;
		}
	}
	
	setSpeed( speed ) {
		this.gridUniforms.waveSpeed.value = speed * 1500.0;
	}
	
	addTime( delta_time ) {
		this.gridUniforms.uTime.value += delta_time * this.gridUniforms.waveSpeed.value;
	}
	
	setTime( tick_time ) {
		this.gridUniforms.uTime.value = tick_time;
	}
}
