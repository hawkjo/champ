r = 5.4;
iw = 136;
ih = 95;
id = 8;

ow = 159.5+1.25;
oh = 109.5+1.25;
od = 3;

sw = 51;
sh = 28;
sd = 6;
st = 4;

hw = 5;
hh = sw+15;
hd = 3;

pw = sw/4;
ph = sh/2;
pd = 3;

//Plate for Tubing Pressure
//translate([70.875, 26.7, -100])
//difference () {
//cube([16.725,42.85,0.5]);
//    
//    translate([9.44,2.5,0])
//    cube([4.15,5,4.2]);
//
//    translate([9.44,35.35,0])
//    cube([4.15,5,4.2]);    
//}

//This code builds the Miseq Microscope Adapter.
union () {
difference () {
		union () {
			union() {
				translate( [0,r,0] ) cube( [iw, ih-2*r, id] );
				translate( [r,0,0] ) cube( [iw-2*r, ih, id] );
				translate( [r,r,0] ) cylinder( r=r, h=id, $fn=32 );
				translate( [iw-r,r,0] ) cylinder( r=r, h=id, $fn=32 );
				translate( [r,ih-r,0] ) cylinder( r=r, h=id, $fn=32 );
				translate( [iw-r,ih-r,0] ) cylinder( r=r, h=id, $fn=32 );
			}
//This code builds the top layer of the adapter.			
			translate( [(iw-ow)/2,(ih-oh)/2,id-od] )
			union() {
				translate( [0,r,0] ) cube( [ow, oh-2*r, od] );
				translate( [r,0,0] ) cube( [ow-2*r, oh, od] );
				translate( [r,r,0] ) cylinder( r=r, h=od, $fn=32 );
 				translate( [ow-r,r,0] ) cylinder( r=r, h=od, $fn=32 );
				translate( [r,oh-r,0] ) cylinder( r=r, h=od, $fn=32 );
				translate( [ow-r,oh-r,0] ) cylinder( r=r, h=od, $fn=32 );
			}

}
//This code builds the screw openings along the adapter.
for( x=[7:8:iw-7] ) {
			translate([x,ih-r-79,-1]) cylinder( d=4, h=15, $fn=30 );
			translate([x,ih-16+5.4,-1]) cylinder( d=4, h=15, $fn=30 );
			translate([x,ih-r-79,-1]) cylinder( d=5.9, h=4, $fn=30 );
			translate([x,ih-16+5.4,-1]) cylinder( d=5.9, h=4, $fn=30 );
}

//Chip Holder	
		translate([-15,((ih-sh)/2),1])
		cube([154,sh,id+od+2]);

		translate([((iw-sw)/2)+42,((ih-sh)/2),-1])
		cube([53.5,28,3]);	

		translate([0,((ih-sh)/2),0])
		cube([28.85+((iw-sw)/2)-24.5,sh,id+od+2]);

	}

}

