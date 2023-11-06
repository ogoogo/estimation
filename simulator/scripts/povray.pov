#include "textures.inc"
// #include "Woods.inc"
#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "skies.inc"
#fopen information "inforow.txt" read
#fopen moon_dcm "moondcm.txt" read

#read(information, id, et, r_x, r_y, r_z, moon_x, moon_y, moon_z, sun_x, sun_y, sun_z, c_e11, c_e12, c_e13, c_e21, c_e22, c_e23, c_e31, c_e32, c_e33, c11, c12, c13, c21, c22, c23, c31, c32, c33, sun_d_x, sun_d_y, sun_d_z,mn11,mn12,mn13,mn21,mn22,mn23,mn31,mn32,mn33 )
// #read(moon_dcm, mn11,mn12,mn13,mn21,mn22,mn23,mn31,mn32,mn33)

#declare L_moon = <moon_x,moon_y,moon_z>;
#declare L_sun = <sun_x, sun_y, sun_z>;

#declare L_equ = <r_x, r_y, r_z>;

global_settings{
    ambient_light <0, 0, 0> // <red, green, blue}
	hf_gray_16 on
}

camera{ 
	// location L_equ
	location <0,0,0>
	sky <0,0,1>
	look_at <0,1,0>

	matrix <c11, c12, c13,
			c21, c22, c23,
			c31, c32, c33, 
			r_x, r_y, r_z>
	// matrix <1,0,0,
	// 		3, 1, 1,
	// 		0, 0, 1, 
	// 		0, 0, 0>
	
	// look_at L_equ + <3, 1, 1>

	// look_at L_equ + Att_equ
	// look_at <0,0,0>
	// look_at L_moon

   
    right -x*image_width/image_height
	angle 5.58

} 
light_source { 
	L_sun
	color rgb <2,2,2>
    parallel
} 


// sky_sphere{
// 	S_Cloud2
// }
object { 
	sphere {<0 , 0 , 0> , 1737.4}
	texture { 
		pigment { 
			// rgb <1 , 1 , 1>
			image_map {
 				jpeg "moonmap.jpg" map_type 1
			} 
		}
		finish{
			diffuse 1.5
			brilliance 2
			crand 0
			// phong 0.6
		}
	} 
	matrix
		<mn11,mn12,mn13,
		 mn21,mn22,mn23,
		 mn31,mn32,mn33,
		 0,0,0>
	

    translate L_moon // <x, y, z>

}
object { 
	sphere {<0 , 0 , 0> , 6378.1} 
	
	texture { 
		pigment { rgb <1 , 1 , 1> } 
	} 

} 

// merge {
// 	object { Cylinder_X scale 1000 pigment { color Red } }		// x軸
// 	object { Cylinder_Y scale 1000 pigment { color Green } }	// y軸
// 	object { Cylinder_Z scale 1000 pigment { color Blue } }		// z軸
// }

// object{
// 	sphere{<100000,0,0>,5000 pigment{color Red}}
	
// }
// object{
// 	sphere{<0,100000,0>,5000 pigment{color Green}}
	
// }
// object{
// 	sphere{<0,0,100000>,5000 pigment{color Blue}}
// }

// object {
//   plane { z, 100000 }
//   pigment { 
//     checker
//       color rgb 1,
//       color rgb <1,0,0>
//   }
// }

// object {
//   plane { z, 100001 }
//   pigment { 
//     checker
//       color rgb 1,
//       color rgb <0,0,1>
//   }
// }
