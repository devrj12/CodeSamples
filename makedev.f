      program readpim
c
      dimension alts(101)
      dimension dens(101)
      dimension dev(201,201,101)
c
      character*80 line
c
      open(1,file="PIM3D.OUT",form="formatted",status="old")
      open(2,file="PIMMAT.dat",form="formatted",status="unknown")
      open(3,file="ALTITUDE.dat",form="formatted",status="unknown")
c
      read(1,100)line
      read(1,100)line
c
  100 format(a80)
c
      read(1,*)iy,id,ut,f107,pkp,ssn
c
      hr=(ut/3600.0)
c
      print101,iy,id,ut,f107,pkp
  101 format("PIM (in MHz) for ",i4,"/",i3,1x,f7.1," UT",
     x       " F10.7 ",f6.1," Kp ",f4.1)
c
      do 10 i=1,5
      read(1,100)line
      print*,line
   10 continue
c
      read(1,*)starter,ender,stlon,enlon,nlats,nlons,dlat,dlon
c
      nalts=101
c      
      do 12 i=1,3
      read(1,100)line
      print*,line
   12 continue
c
      read(1,*)alts
c
      do 13 i=1,101
   13 write(3,*)alts(i)
c
c   read in the densities
c
      print*,"Lats:",nlats
      print*,"Lons:",nlons
      print*,"Alts:",nalts
c
      do 14 ilat=1,nlats
      do 14 ilon=1,nlons
c
      read(1,100)line
      read(1,100)line
      read(1,100)line
c
      read(1,*)dens
c
      do 15 ialt=1,nalts
   15 dev(ilat,ilon,ialt)=dens(ialt)
c
   14 continue
c
      do 20 ilat=1,nlats
      do 20 ilon=1,nlons
      do 20 ialt=1,nalts
c
      write(2,201)dev(ilat,ilon,ialt)
  201 format(1x,1pe10.2)
c
   20 continue
c
      stop
      end
