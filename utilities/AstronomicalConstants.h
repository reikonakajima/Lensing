// 	$Id: AstronomicalConstants.h,v 1.7 2006-07-07 18:47:51 garyb Exp $
// Include file with common "magic numbers" for Astronomy
#ifndef ASTROCONST_H
#define ASTROCONST_H

#ifndef PI
#define PI 3.14159265358979323
#endif

// conversion constants:
const double RadToArcsec=3600.*180./PI;
const double AU=1.49597870691e11;	   //Astronomical Unit, m
const double Parsec=AU*RadToArcsec;
const double ABZeropoint=3.63e-23;         //AB=0 mag, in W/Hz/m^2
const double Jansky=1e-26;	           //convert Jy to W/Hz/s/m^2
const double Planck=6.626068e-34;          //Planck's const in Js.
const double SpeedOfLight=2.99792458e8;    //speed of light in m/s
const double Boltzmann=1.38065e-23;	   //Boltzmann's const, J/K
const double Micron=1e-6;		   //micron, in m
const double HubbleLengthMpc=SpeedOfLight*1e-5; //c/H0, in h^-1 Mpc
const double HubbleDistance=HubbleLengthMpc*1e6*Parsec; //c/H0, in h^-1 m
const double RecombinationRedshift=1088.;

// Orbital Mechanics constants, units of AU, solar mass, and Year
const double  TPI           = 2.*PI;
const double  DEGREE	    = PI/180.;               // One degree, in rad
const double  GM            = 4.*PI*PI/1.0000378;    //solar gravitation
const double  SolarSystemMass= 1.00134;
const double  ARCSEC	    = PI/180./3600.;
const double  ARCMIN	    = PI/180./60.;
const double  YEAR          = 1.;
const double  DAY	    = 1./365.25;	     //Julian day, 86400 s
const double  HOUR          = DAY/24.;
const double  MINUTE        = HOUR/60.;
const double  SECOND        = MINUTE/60.;
const double  SpeedOfLightAU= 63241.06515;	     //in AU/YR

//Obliquity of ecliptic at J2000
const double  EclipticInclination =23.43928*DEGREE;  
const double  EclipticNode        =0.;
//Inclination and ascending node of invariable plane in the J2000
//equatorial system: ??? check this
const double  InvariableInclination=(23+22.11/3600.)*DEGREE;
const double  InvariableNode       = (3+(52.+23.7/60.)/60.)*DEGREE;

const double EarthMass        = 3.00349e-6;	//Earth mass in Solar units
const double MJD0             =2400000.5;	//Offset for modified Julian dates


#endif  // ASTROCONST_H

