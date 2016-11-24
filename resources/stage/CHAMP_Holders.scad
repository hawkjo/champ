//MiSeq chip
module MiSeqChip(){
color("red")
    translate([0,0,-1])
difference(){
    cube([28,14,1.26]);
    translate([1,7,0])
        cylinder(r=0.25,h=1,$fn=32);
    translate([1,5,0])
        cylinder(r=0.25,h=1,$fn=32);
    translate([0.75,2.5,1])
        cube([21,7,1.26]);
}}

//MiSeq gasket
module MiSeqGasket(){
    color("blue")
difference(){
union(){
    cube([17,21,2.1]);
    translate([0,2.5,0])
        cube([23.5,16,2.1]);
}
translate([21,3.5,1])
    cube([2.5,14,1.1]);
}}

module MiSeqPressurePlate(){
    translate([0,30,0])
union(){
    difference() {
        translate([-42,1.5,0])        //base
            cube([84,38.5,4]);
        //screw openings
        translate([36.2,6,0])
            cylinder(r=2.5,h=6,$fn=16);
        translate([-36.2,6,0])
            cylinder(r=2.5,h=6,$fn=16);
        translate([36.2,30,0])
            cylinder(r=2.5,h=6,$fn=16);
        translate([-36.2,30,0])
            cylinder(r=2.5,h=6,$fn=16);
        //interface opening
        translate([-12,0,0]) 
            cube([24,29,6]);
    //fitting holes 
        translate([-7,35,0])         
            cylinder(r=3.6,h=6,$fn=6);
        translate([7,35,0])
            cylinder(r=3.6,h=6,$fn=6);
    }
    //gasket holder  
    difference(){
        //base  
        translate([-13.8,4,-10.5])  
            cube([27.6,13.5,14.5]);
        //gasket slot
        translate([-10.5,6.5,-10.5])
            cube([21,17.5,1.9]);
        translate([-8,1.5,-10.5])
            cube([16,17.5,1.9]);

        }
}}

//MiSeqCradle
module MiSeqCradle (){
color ("purple")
    translate([0,30,0])
union() {
    difference() {
        union(){
            translate([-42,0,8])      //base left
                cube([38,44,3]);
            translate([4,0,8])        //base right
                cube([38,44,3]);
            }
        translate([36.2,10,7])        //screw openings
            cylinder(r=2.5,h=6,$fn=16);
        translate([-36.2,10,7])
            cylinder(r=2.5,h=6,$fn=16);
        translate([36.2,34,7])
            cylinder(r=2.5,h=6,$fn=16);
        translate([-36.2,34,7])
            cylinder(r=2.5,h=6,$fn=16);
        translate([-6.2,0,6])       //prism opening
            cube([12.4,32,9]);
        translate([-15,12,7])       //heater screw openings
            cylinder(r=2.2,h=6,$fn=16);
        translate([15,12,7])       
            cylinder(r=2.2,h=6,$fn=16);   
    }
    difference() {
        translate([-12.2,0,-3])
            cube([24.4,44,12]);
        translate([-6,0,-4])        //scope opening
            cube([12,25,12]);
        translate([0,25,-4])
            cylinder(r=6,h=6,$fn=16);
        translate([-7,0,-2]) 
            cube([14,25,1.3]);      //chip notch
        translate([-6.2,0,-0.7]) 
            cube([12.4,32,9.7]);    //prism opening
        translate([-4,28,-0.9])     //top laser bevel
            rotate([30,0,0]){
                cube([8,35,7.7]);}
        }  
}}

translate([23.5,10.5,9])
    rotate([0,0,90]){
    MiSeqPressurePlate();}
translate([21,3.5,0.5])
    MiSeqChip();
translate([0,0,-1.5])
    MiSeqGasket();
translate([24,10.5,2.5])
    rotate([0,0,-90]){
        MiSeqCradle();}


