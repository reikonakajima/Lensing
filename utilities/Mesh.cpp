#include "Mesh.h"


/*
template <class Ttype, class Tpos> 
Mesh<Ttype, Tpos>::Mesh(int nn1, int nn2, int nn3,
			vector<Ttype>& P, bool period) {
  isPeriodic=period;
  nm.resize(3);
  nm[0] = nn1;
  nm[1] = nn2;
  nm[2] = nn3;
  int np= P.size();
  head.resize(nn1*nn2*nn3);
  next.resize(np);
  vector<int> tail(nn1*nn2*nn3);
  dat = &P;
  for (int nn=0; nn<nn1*nn2*nn3; nn++) head[nn]=tail[nn]=-1;
  for (int nn=0; nn<np; nn++) next[nn]=-1;
  for (int nn=0; nn<np; nn++) {
    int ix = int(nn1*periodic(P[nn].pos[0]));
    int iy = int(nn2*periodic(P[nn].pos[1]));
    int iz = int(nn3*periodic(P[nn].pos[2]));
    int ii = nn2*nn3*ix+nn3*iy+iz;
    if (head[ii]<0) {	// Cell is empty.
      head[ii]=tail[ii]=nn;
    }
    else {
      next[tail[ii]] = nn;
      tail[ii] = nn;
    }
  }
}
*/


template <class Ttype, class Tpos> 
Mesh<Ttype, Tpos>::Mesh(int nn1, int nn2, int nn3,
			vector<Ttype>& P, bool period,
			Tpos xmin_, Tpos xmax_,
			Tpos ymin_, Tpos ymax_,
			Tpos zmin_, Tpos zmax_) : 
  xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), zmin(zmin_), zmax(zmax_) 
{
  isPeriodic=period;
  nm.resize(3);
  nm[0] = nn1;
  nm[1] = nn2;
  nm[2] = nn3;
  dx = (xmax-xmin) / nn1;
  dy = (ymax-ymin) / nn2;
  dz = (zmax-zmin) / nn3;
  xmid = (xmax+xmin) / 2.0;
  ymid = (ymax+ymin) / 2.0; 
  zmid = (zmax+zmin) / 2.0;
  int np= P.size();
  head.resize(nn1*nn2*nn3);
  next.resize(np);
  vector<int> tail(nn1*nn2*nn3);
  dat = &P;
  for (int nn=0; nn<nn1*nn2*nn3; nn++) head[nn]=tail[nn]=-1;
  for (int nn=0; nn<np; nn++) next[nn]=-1;
  for (int nn=0; nn<np; nn++) {
    int ix = int((periodicX(P[nn]->getX())-xmin)/dx);
    int iy = int((periodicY(P[nn]->getY())-ymin)/dy);
    int iz = int((periodicZ(P[nn]->getZ())-zmin)/dz);
    int ii = nn2*nn3*ix+nn3*iy+iz;
    if (head[ii]<0) {	// Cell is empty.
      head[ii]=tail[ii]=nn;
    }
    else {
      next[tail[ii]] = nn;
      tail[ii] = nn;
    }
  }
}


template <class Ttype, class Tpos> 
list<int> 
Mesh<Ttype, Tpos>::getNearMeshList(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin) {
  std::list<int> nbr;   // the return list
  nbr.clear();          // clear list
  int ix,iy,iz,ii,p;
  Tpos dst2;
  ix = int((x - xmin)/dx);     // calculate index of point
  iy = int((y - ymin)/dy);              
  iz = int((z - zmin)/dz);
  int srx = int(rmax/dx)+1;    // calculate radius in terms of mesh indicies
  int sry = int(rmax/dy)+1;
  int srz = int(rmax/dz)+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (int i = 0; i < close.size(); ++i)
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {
	dst2 = distance2((*dat)[p]->getX(), (*dat)[p]->getY(), (*dat)[p]->getZ(), x, y, z);
	if (dst2<rmax*rmax && dst2>=rmin*rmin) {
	  nbr.push_back(p);
	}
      } while( (p=next[p])>=0 );
    }
  }
  return(nbr);
}


template <class Ttype, class Tpos> 
multimap<double, int> 
Mesh<Ttype, Tpos>::getNearMeshMap(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin) {
  std::multimap<double, int> nbr;   // the return list
  nbr.clear();          // clear list
  int ix,iy,iz,ii,p;
  Tpos dst2;
  ix = int((x - xmin)/dx);     // calculate index of point
  iy = int((y - ymin)/dy);              
  iz = int((z - zmin)/dz);
  int srx = int(rmax/dx)+1;    // calculate radius in terms of mesh indicies
  int sry = int(rmax/dy)+1;
  int srz = int(rmax/dz)+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (int i = 0; i < close.size(); ++i) {
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {
	dst2 = distance2((*dat)[p]->getX(), (*dat)[p]->getY(), (*dat)[p]->getZ(), x, y, z);
	if (dst2<rmax*rmax && dst2>=rmin*rmin) {
	  nbr.insert(std::pair<double,int>(dst2,p));
	}
      } while( (p=next[p])>=0 );
    }
  }
  }
  return(nbr);
}


template <class Ttype, class Tpos> 
multimap<double, int> 
Mesh<Ttype, Tpos>::getNearAngleMap(Tpos ra, Tpos dec, Tpos z, Tpos thetamax, Tpos thetamin) {
  std::multimap<double, int> nbr;    // the return list
  nbr.clear();                       // clear list
  int ix,iy,iz,ii,p;
  Tpos angsep, cosangsep;
  Tpos cosa, sina, cosb, sinb, cosC;
  ix = int((ra - xmin)/dx);          // calculate index of point
  iy = int((dec - ymin)/dy);              
  iz = int((z - zmin)/dz);
  Tpos cosd_1 = cos((dec+thetamax)*DEGREE);
  Tpos cosd_2 = cos((dec-thetamax)*DEGREE);
  Tpos cosd = cosd_1;
  if ((cosd_1 > cosd_2)) 
    cosd = cosd_2;
  //cerr << "! " << cosd << " " << cosd_1 << " " << cosd_2 << endl;
  //cerr << thetamax << " (" << ra << " " << dec << ")" << endl;
  // XXX BUG NOTE: THIS IS PROBLEMATIC WHEN (dec +/- thetamax) IS > 90 deg !!!  XXXX
  int srx = int(thetamax/dx/cosd)+1; // calculate radius in terms of mesh indicies
  int sry = int(thetamax/dy)+1;      // (ra shrinks by cos(dec), include this effect)
  int srz = int(thetamax/dz)+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {                           // calculate spherical angle separation
	cosb = cos(dec*DEGREE);
	sinb = sin(dec*DEGREE);
	cosa = cos((*dat)[p]->getY()*DEGREE);
	sina = sin((*dat)[p]->getY()*DEGREE);
	cosC = cos((ra - (*dat)[p]->getX()) * DEGREE);
	cosangsep = sina*sinb + cosa*cosb*cosC;
	//cerr << cosangsep << "=" << sina << "*" << sinb << "+" 
	//     << cosa << "*" << cosb << "*" << cosC << endl;
	angsep = acos(cosangsep) / DEGREE;
	//cerr << angsep * 3600 << endl;  // DEBUG
	//return(nbr);  // DEBUG
	if (angsep < thetamax && angsep >= thetamin) {
	  nbr.insert(std::pair<double,int>(angsep,p));
	}
      } while( (p=next[p])>=0 );
    }
  }
  return(nbr);
}


/*
template <class Ttype, class Tpos> 
list<int> 
Mesh<Ttype, Tpos>::nearobj(Tpos pos[], Tpos rmin, Tpos rmax) {
  std::list<int> nbr;
  nbr.clear();
  int ix,iy,iz,ii,p;
  Tpos dst2;
  ix = int(nm[0]*pos[0]);
  iy = int(nm[1]*pos[1]);
  iz = int(nm[2]*pos[2]);
  int srx = int(rmax*nm[0])+1;
  int sry = int(rmax*nm[1])+1;
  int srz = int(rmax*nm[2])+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {
	dst2 = distance2((*dat)[p].pos,pos);
	if (dst2<rmax*rmax && dst2>rmin*rmin) {
	  nbr.push_back(p);
	}
      } while( (p=next[p])>=0 );
    }
  }
  return(nbr);
}


template <class Ttype, class Tpos> 
list<Tpos> 
Mesh<Ttype, Tpos>::dist2(Tpos pos[], Tpos rmin, Tpos rmax) {
  // Returns a list of distances^2 for objects near pos.
  std::list<Tpos> nbr;
  nbr.clear();
  int ix,iy,iz,ii,p;
  Tpos dr,dst2;
  ix = int(nm[0]*pos[0]);
  iy = int(nm[1]*pos[1]);
  iz = int(nm[2]*pos[2]);
  int srx = int(rmax*nm[0])+1;
  int sry = int(rmax*nm[1])+1;
  int srz = int(rmax*nm[2])+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {
	dst2=0;
	dst2 = distance2((*dat)[p].pos,pos);
	if (dst2<rmax*rmax && dst2>rmin*rmin) {
	  nbr.push_back(dst2);
	}
      } while( (p=next[p])>=0 );
    }
  }
  return(nbr);
}
*/


template <class Ttype, class Tpos> 
list<Tpos> 
Mesh<Ttype, Tpos>::getDistSqList(Tpos x, Tpos y, Tpos z, Tpos rmax, Tpos rmin) {
  // Returns a list of distances^2 for objects near pos.
  std::list<Tpos> nbr;
  nbr.clear();
  int ix,iy,iz,ii,p;
  Tpos dr,dst2;
  ix = int((x - xmin)/dx);     // calculate index of point
  iy = int((y - ymin)/dy);              
  iz = int((z - zmin)/dz);
  int srx = int(rmax/dx)+1;    // calculate radius in terms of mesh indicies
  int sry = int(rmax/dy)+1;
  int srz = int(rmax/dz)+1;
  std::vector<int> close = closemeshes(ix,iy,iz,srx,sry,srz);
  for (std::vector<int>::iterator ii=close.begin(); ii!=close.end(); ii++) {
    if ( (p=head[*ii])>=0 ) {
      do {
	dst2=0;
	dst2 = distance2((*dat)[p]->getX(), (*dat)[p]->getY(), (*dat)[p]->getZ(), x, y, z);
	if (dst2<rmax*rmax && dst2>rmin*rmin) {
	  nbr.push_back(dst2);
	}
      } while( (p=next[p])>=0 );
    }
  }
  return(nbr);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::periodicX(Tpos x) {	// Wraps x periodically in the range [min,max).
  Tpos tmp = x;
  while (tmp >= xmax) tmp = x-dx;
  while (tmp <  xmin) tmp = x+dx;
  return(tmp);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::periodicY(Tpos y) {	// Wraps x periodically in the range [min,max).
  Tpos tmp = y;
  if (tmp >= ymax) tmp = y-dy;
  if (tmp <  ymin) tmp = y+dy;
  return(tmp);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::periodicZ(Tpos z) {	// Wraps x periodically in the range [min,max).
  Tpos tmp = z;
  if (tmp >= zmax) tmp = z-dz;
  if (tmp <  zmin) tmp = z+dz;
  return(tmp);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::diffperiodicX(Tpos x) {	// Periodic differences.
  Tpos   tmp;
  tmp = x;
  do {
    if (tmp>=xmid) tmp -= dx;
    if (tmp<-xmid) tmp += dx;
  } while( tmp>=xmid || tmp<-xmid);
  return(tmp);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::diffperiodicY(Tpos y) {	// Periodic differences.
  Tpos   tmp;
  tmp = y;
  do {
    if (tmp>=ymid) tmp -= dy;
    if (tmp<-ymid) tmp += dy;
  } while( tmp>=ymid || tmp<-ymid);
  return(tmp);
}


template <class Ttype, class Tpos> 
Tpos 
Mesh<Ttype, Tpos>::diffperiodicZ(Tpos z) {	// Periodic differences.
  Tpos   tmp;
  tmp = z;
  do {
    if (tmp>=zmid) tmp -= dz;
    if (tmp<-zmid) tmp += dz;
  } while( tmp>=zmid || tmp<-zmid);
  return(tmp);
}



/*
template <class Ttype, class Tpos> 
Tpos	
Mesh<Ttype, Tpos>::distance2(Tpos a[], Tpos b[]) {
  // Returns the distance^2 between points a and b.
  Tpos dr,dst2=0;
  for (int j=0; j<3; j++) {
    if (isPeriodic==true)
      dr  = diffperiodic(a[j]-b[j]);
    else
      dr  = a[j]-b[j];
    dst2 += dr*dr;
  }
  return(dst2);
}
*/



template <class Ttype, class Tpos> 
Tpos	
Mesh<Ttype, Tpos>::distance2(Tpos ax, Tpos ay, Tpos az, Tpos bx, Tpos by, Tpos bz) {
  // Returns the distance^2 between points a and b.
  Tpos dx, dy, dz, dst2=0;
  if (isPeriodic==true) {
    dx  = diffperiodicX(ax-bx);
    dy  = diffperiodicY(ay-by);
    dz  = diffperiodicZ(az-bz);
  } else {
    dx  = ax-bx;
    dy  = ay-by;
    dz  = az-bz;
  }
  dst2 = dx*dx + dy*dy + dz*dz;
  return(dst2);
}



template <class Ttype, class Tpos> 
vector<int> 
Mesh<Ttype, Tpos>::closemeshes(int ix, int iy, int iz,
			       int srx, int sry, int srz) {
 // Returns a list (vector) of the mesh points within sr of i.
  int ii,nn;
  vector<int> retlist;
  retlist.clear();
  nn = 0;
  if (isPeriodic==true) {
    retlist.resize( (2*srx+1)*(2*sry+1)*(2*sry+1) );
    for (int iix=-srx; iix<=srx; iix++)
      for (int iiy=-sry; iiy<=sry; iiy++)
	for (int iiz=-srz; iiz<=srz; iiz++) {
	  ii = nm[1]*nm[2]*( (ix+iix+nm[0])%nm[0] ) +
	    nm[2]*( (iy+iiy+nm[1])%nm[1] ) +
	    ( (iz+iiz+nm[2])%nm[2] ) ;
	  retlist[nn] = ii;
	  nn++;
	}
  }
  else {
    //cerr << srx << endl;
    int srxmin = ix-srx;  if (srxmin<  0   ) srxmin=0;  if (srxmin>=nm[0]) srxmin=nm[0]-1;
    //cerr << srxmin << endl;
    int srymin = iy-sry;  if (srymin<  0   ) srymin=0;  if (srymin>=nm[1]) srymin=nm[1]-1;
    //cerr << srymin << endl;
    int srzmin = iz-srz;  if (srzmin<  0   ) srzmin=0;  if (srzmin>=nm[2]) srzmin=nm[2]-1;
    //cerr << srzmin << endl;
    int srxmax = ix+srx;  if (srxmax>=nm[0]) srxmax=nm[0]-1;  if (srxmax<  0   ) srxmax=0; 
    //cerr << srxmax << endl;
    int srymax = iy+sry;  if (srymax>=nm[1]) srymax=nm[1]-1;  if (srymax<  0   ) srymax=0; 
    //cerr << srymax << endl;
    int srzmax = iz+srz;  if (srzmax>=nm[2]) srzmax=nm[2]-1;  if (srzmax<  0   ) srzmax=0; 
    //cerr << srzmax << endl;
    retlist.resize( (srxmax-srxmin+1)*(srymax-srymin+1)*(srzmax-srzmin+1) );
    //cerr << retlist.size() << endl;
    nn = 0;
    for (int iix=srxmin; iix<=srxmax; iix++)
      for (int iiy=srymin; iiy<=srymax; iiy++)
	for (int iiz=srzmin; iiz<=srzmax; iiz++) {
	  ii = nm[1]*nm[2]*iix+nm[2]*iiy+iiz;
	  retlist[nn] = ii;
	  nn++;
	}
    //cerr << nn << endl;
  }
  return(retlist);
}

#include "LensObjects.h"
template class Mesh<LensObject*, double>;    // explicit instantiation
