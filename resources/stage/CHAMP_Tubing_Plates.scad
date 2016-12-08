//translate([49,25,8])
union () {
difference () {
cube([33.45,42.85,1.5]);

//Screw Cubes
//top left
translate([3.1,-0.1,-0.1])
cube([4.15,7.5,4.2]);

//bottom left
translate([26.01,-0.1,-0.1])
cube([4.15,7.5,4.2]);

//top right
translate([3.1,35.35,-0.1])
cube([4.15,7.5,4.2]);

//bottom right
translate([26.01,35.35,-0.1])
cube([4.15,7.5,4.2]);

rotate([0,90,0])
translate([-1.5,21.75,2])
cylinder(d=1.25, h=34, $fn=32);

rotate([0,90,0])
translate([-1.5,23.25,2])
cylinder(d=1.25, h=34, $fn=32);

translate([2,21.75,-1.3])
cylinder(d=1.25, h=4.2, $fn=32);

translate([2,23.25,-1.3])
cylinder(d=1.25, h=4.2, $fn=32);
}

translate([0.625,20.375,-3.4])
difference () {
cube([2.5,4,4]);

translate([1.375,1.375,-.1])
cylinder(d=1.25, h=4.2, $fn=32);

translate([1.375,2.875,-.1])
cylinder(d=1.25, h=4.2, $fn=32);


}
}


//Plate for Tubing Pressure
translate([16.725,0, -4])
difference () {
cube([16.725,42.85,0.5]);
    
    translate([9.44,2.5,0])
    cube([4.15,5,4.2]);

    translate([9.44,35.35,0])
    cube([4.15,5,4.2]);    
}
