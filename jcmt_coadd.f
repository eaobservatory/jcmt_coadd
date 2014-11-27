      PROGRAM   JCMT_COADD

c History:
c   written: David Hughes (Lancs. Poly, LPVAD::DHH) 01/06/91
c
c   05-Aug-1991: (JACH::FJO) get data directly from GSD data file 
c   05-Aug-1991: (JACH::FJO) compile and link with following:
c   25-Oct-1991: (JACH::FJO) handle Cycle Reversal TRUE or FALSE
c   29-Dec-1991: (LPVAD::DHH)
c   08-Jan-1991: (LPVAD::DHH) corrects for variable extinction, modified
c                plotting routine, allows user to repeat extinction correction,
c                despiking or K-S test, calibrates data if sensitivity known,
c                writes summary of reduction to output file.
c   05-Mar-1992: (JACH::FJO) Check ERRDATA(i) for spikes and fix DATA(i).
c   09-Mar-1992: (JCMT::FJO) Automatic terminal switching from Tek to VT300
c   30-Apr-1992: (JACH::FJO) Combine files with even or odd cycles,
c                            cycle reversal or no cycle reversal.
c   24-Feb-1993: (JCMT::FJO) Increase amstart, amend, indiviobs, and keepscan
c                            from 20 to 100 (allows 100 obs to be co-added)
c   05-Mar-1994: (JACH::RPT) Line parser, repeat option 4, hardcopy option.
c                            CBE files and inquire (HEM)
c   30-Jun-1994: (JACH::FJO) Read aperture as real or integer

C ------------- UNIX version -----------------

c   10-Jan-1996: (JACH:TIMJ) Port to unix

C   Changes for UNIX:

C     0) Wrote simple makefile
C     1) GSD problems:
C               Filename should not have .dat in it
C               Data TYPEs are now upper case NOT lower case
C               Use GSD_INQ_SIZE with GSD__MXDIMS dimensions
C     2) Files:
C               Shouldn't try to read from unit 6!
C               Lots of errors caused by OPEN command specyfying 'NEW'
C     3) GET_CHARACTER modified to be standalone. Had trouble passing
C        Strings back through (real)rvalues array. Gave up and did
C        own gsd_get0c (meant that opening and closing GSD file
C        is necessary 3 times more per file :-( - don't have time to
C        be tidy!)

C  13 Feb 1997
C     Fix filename inquire for CBE data
C     Fix internal read for aperture if no aperture (ie CBE)
C     Remove NAG (replace KS test with SCUBA/Num Rec version)

C ------------ End UNIX fixes ----------------

*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 3 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with this program.  If not, see <http://www.gnu.org/licenses/>.

*  Copyright
*     1991-1997 Science and Engineering Research Council.
*     All Rights Reserved.

c End History
c
c  $!COMPILE_COADD.COM
c  $ IF F$TRNLNM("MT_ROOT") .EQS. "" THEN $JCMTSETUP
c  $ FORTRAN/EXTEND/OPTIMIZE JCMT_COADD
c  $ EXIT
c  $!LINK_COADD.COM
c  $ IF F$TRNLNM("MT_ROOT") .EQS. "" THEN $JCMTSETUP
c  $ define/user ems_dir lstardisk:[starlink.lib.ems.standalone]
c  $ define/user pgplot_dir lstardisk:[starlink.lib.pgplot]
c  $ define/user gns_dir lstardisk:[starlink.lib.gns]
c  $ define/user chr_dir lstardisk:[starlink.lib.chr]
c  $ define/user psx_dir lstardisk:[starlink.lib.psx]
c  $ define/user cnf_dir lstardisk:[starlink.lib.cnf]
c  $ define/user gks_dir lstardisk:[starlink.lib.gks]
c  $ LINK JCMT_COADD,-
c    nag_lib/lib,-
c    pgplot_dir:grpckg/lib,-
c    gns_dir:gns/lib,-
c    chr_dir:chr/lib,-
c    psx_dir:psx_link/opt,-
c    cnf_dir:cnf_link/opt,-
c    mt_gsddir:gsd_lib/opt,-
c    gks_dir:gkslink/opt,-
c    @ems_dir:emslink
c  $ EXIT

      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'

*    Local variables :
	REAL MVOLTS_ZERO(3000),MVOLTS(3000), PAIR(3000)
        REAL MVLEFT(3000),MVRIGHT(3000), ERR_ZERO(3000)
        REAL ERRRIGHT(3000), ERRLEFT(3000), ERR(3000)
	REAL MEANFIN, MEANMV(3000), FIRSTERROR, STDEV
        REAL ERRORFIN, ERRORSTDEV, MEANERROR(3000)
        REAL A(3000), B(3000), OBS(3000)
        REAL MV_MEDIAN, aperture, signoise
        character*8 filter
        real ra,dec,lat,lst,el,lststart,lstend
        real z1,z2,h,airmass1,airmass2,trig
        real amstart(100),amend(100), overheads

        INTEGER DELPAIR(300), NODEL, DATAFILE, iaperture
	INTEGER NOBS, N, M, COUNTER, NEWOBS, DUMMY, ENDFILE
        integer n1,diff,k,j,i,nitems,isave,start,end
	integer samplesize,totalint,c, indivobs(100), nocoadds

	CHARACTER*20 OUTPUT, sourcename
	character*80 variname, cdummy
        character*3 ans
        character*4 scan
        integer keepscan(100)
        character*10 stringobs

*  New variables: created by FJO 08-Sep-1991 and 25-Oct-1991
        character*16 c_filter, c_aperture
        REAL temp_array(800)
        REAL*8 data(400), errdata(400)
        REAL*8 dvalue
        integer psign
        integer n_data
        integer STATUS
        character*64 input, datadir
        character*128 actinput
        logical cycle_reversal, end_even, new_even

*  New variables: created by DHH 04-Jan-1992
	real taumean, taustdev, respons, calmean, calerror, stderror
	integer savenobs, ksnobs, option
 	character*4 dans, kans, optionchar
	character*20 summary

*  New variables: added by RPT 05-Mar-1994
        integer      cline, length
        logical      exists,nofile     ! HEM

* Check length of string with system library call: TimJ 13-Feb-1997
        integer lnblnk
        external lnblnk

1000    FORMAT (A)
1001    format(a8)     !  filter        
1002    format(i6)     !  ncycles
1003    format(a20)    !  source
1004    format(f5.1)   !  aperture
1010	format(4x,'Source',2x,a20,2x,'Filter',2x,a8,
     :          2x,'Aperture',2x,f5.1)
2001  format(f23.20)
2002  format(i6)

C  Get DATADIR from the command line parameters
C 
        call getenv( "DATADIR", datadir )
        length = 64
        do while ( datadir(length:length).eq.' ' .and. length.gt.0 )
          length = length -1
        enddo

	write(6,*)' '
	write(6,*)' '
	write(6,'(A)')' Summary of reduction recorded in file '
	write(6,'(A)')'$                      [ summary.dat ] :  '
	read(5,'(A)') summary
	if (summary(1:1) .eq. ' ') then
	  summary = 'summary'
	  open(unit=50, file='summary.dat', status='unknown')	  
	else 
	  open(unit=50, file=summary, status='unknown')	  
	end if


        ans = 'n'
        end = 0
        c = 0
        cline = 0

        write(6,*)' '
        write(6,*)' '
 990    call parseline(keepscan,cline)

 995    c = c + 1                                ! Each scan in turn
        i = keepscan(c)
        nofile= .true.
        scan = '0000'
        if ( i .ge. 1000) then
          write(scan(1:4),'(I4)') i
        else if (i .ge. 100) then
          write(scan(2:4),'(I3)') i
        else if (i .ge. 10) then
          write(scan(3:4),'(I2)') i
        else
          write(scan(4:4),'(I1)') i
        endif

C The GSD library opens a file with or without the .dat extension
C Inquire must have the .dat.

        input='obs_ukt14_'//scan(1:4)//''
C This is a kludge - should use LEN of some kind
        actinput=input(1:14)//'.dat'//''
C End kludge
C        print *,'Filename is ',actinput
        inquire(file=actinput,exist=exists)
        if(exists) then
           nofile=.false.
        else
C Same kludge here...
           actinput=datadir(1:length)//'/'//input(1:14)//'.dat'//''
C        print *,'Filename is ',actinput
           inquire(file=actinput,exist=exists)
           if(exists) then
              nofile=.false.
           else
              input='obs_cbe_'//scan(1:4)//''
C Same kludge here...
              actinput=input(1:12)//'.dat'//''
C        print *,'Filename is ',actinput
              inquire(file=actinput,exist=exists)
              if(exists) then
                 nofile=.false.
                 write(6,*) '   Heterodyne (CBE) photometry:'
                 write(50,*) '   Heterodyne (CBE) photometry:'
              else
                 input='obs_cbe_'//scan(1:4)//''
C Same kludge here...
                 actinput=datadir(1:length)//'/'//input(1:12)//'.dat'//''
C        print *,'Filename is ',actinput
                 inquire(file=actinput,exist=exists)
                 if(exists) then
                    nofile=.false.
                    write(6,*) '   Heterodyne (CBE) photometry:'
                    write(50,*) '   Heterodyne (CBE) photometry:'
                 endif
              endif
           endif
        endif
        if(nofile) then
           write(6,*) '   WARNING: No such UKT or CBE scan ',i
           call compress(keepscan,cline,c)
           goto 998
        endif

        write(6,'(''    File: '',A)') input

	write(50,*) ' '
        write(50,'(3x,a64)') input 
	write(50,*) ' '

           status = adam__ok

           variname='C1SNA1'
           CALL GET_CHARACTER(input, variname, sourcename, STATUS)
           if (STATUS.ne.adam__ok) sourcename=' '
*           write(6,*) 'C1SNA1 =', sourcename

           variname='C4ERA'
           STATUS = ADAM__OK
           CALL GET_DOUBLE(input, variname, dvalue, STATUS) 
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Ra value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
           ra=dvalue
*           write(6,*) 'C4ERA ', ra
  
           variname='C4EDEC'
           STATUS = ADAM__OK
           CALL GET_DOUBLE(input, variname, dvalue, STATUS)
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Dec value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
           dec=dvalue
*           write(6,*) 'C4EDEC ', dec

           variname='C4EL'
           STATUS = ADAM__OK
           CALL GET_DOUBLE(input, variname, dvalue, STATUS) 
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading EL value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
           el=dvalue
*           write(6,*) 'C4EL ', el       

           variname='C1LAT'
           STATUS = ADAM__OK
           CALL GET_DOUBLE(input, variname, dvalue, STATUS)
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Lat value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
           lat=dvalue
*           write(6,*) 'C1LAT ', lat

           variname='C3LST'
           STATUS = ADAM__OK
           CALL GET_DOUBLE(input, variname, dvalue, STATUS)
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading LST value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
           lst=dvalue
*            write(6,*) 'C3LST ', lst

           variname='C3SRT'           
           STATUS = ADAM__OK
           CALL GET_INTEGER(input, variname, totalint, STATUS)
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading SRT value from file'
	      call compress(keepscan,cline,c)
              go to 998
           endif
*            write(6,*) 'C3SRT ', totalint

           variname='C7FIL'
           STATUS = ADAM__OK
           CALL GET_CHARACTER(input, variname, c_filter, STATUS)
*           write(6,*) 'C7FIL =', c_filter
           if (status.ne.adam__ok) c_filter=' '
           filter = c_filter

           variname='C7AP'
           STATUS = ADAM__OK
           CALL GET_CHARACTER(input, variname, c_aperture, STATUS)
*           write(6,*) 'C7AP =', c_aperture, status

* sometimes the filter is in integer form: e.g. 65 instead of 65.0
           if (status.ne.adam__ok) then
              c_aperture='0'
              aperture = 0.0
           else
              if (index(c_aperture,'.') .ne. 0) then
                 read(c_aperture,*) aperture                 
              else
                 if (lnblnk(c_aperture) .gt. 0) then
                    read(c_aperture,*) iaperture
                    aperture=float(iaperture)
                 else
                    aperture = 0.0
                 endif
              endif
           endif

           variname='C3NCYCLE'
           STATUS = ADAM__OK
           CALL GET_INTEGER(input, variname, nitems , STATUS)
           if (STATUS.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Ncycle value from file'
	      call compress(keepscan,cline,c)
              go to 998
           else
*              write(6,*) 'C3NCYCLE ', nitems
           endif

           variname='C6CYCLREV'
           STATUS = ADAM__OK
           CALL GET_LOGICAL(input, variname, cycle_reversal, STATUS)
  
*   set n_data, max is 400 for array d_data
           n_data = 2*nitems
           if (n_data.gt.400) then
              n_data = 400
              nitems = n_data/2
              write(6,*) 'There are more than 400 maximum data values.'
              write(6,*) 'Only the first 400 will be used.'
           endif
           STATUS = ADAM__OK
           CALL GET_C13SPV(input, temp_array, data, n_data, STATUS)
           if (status.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Spectral values from file'
	      call compress(keepscan,cline,c)
              go to 998
           elseif (mod(n_data,2).ne.0) then
              write(6,*) 'Missing Data - ON/OFF mismatch'
              write(6,*) 'Trying to continue'
           endif

*  Fix up nitems if necessary
           if (n_data.lt.2*nitems) then
              nitems = n_data/2
              write(6,*) 'There is less data items than expected'
           write(6,*) 'Trying to continue with ', nitems, ' of paired data'
           endif

           write(6,1010) sourcename, filter, aperture
	   write(50,*) ' '
           write(50,1010) sourcename, filter, aperture
	   write(50,*) ' '

*           write(6,*) 'Spectral Values:'
*           write(6,'(3(e20.10,4x)))') (data(i),i=1,n_data)
*           write(6,*)
*           write(20,*) 'File:', input,
*     :                 ' Source: ', sourcename,
*     :                 'Fil: ', filter,
*     :                 'Ap: ', aperture
*           write(20,*) 'Spectral Values:'
*           write(20,'(3(e20.10,4x)))') (data(i),i=1,n_data)
*           write(20,*)
	        
*   Value of n_data was set above
           STATUS = ADAM__OK
           CALL GET_C13RAW_ERROR(input, temp_array, errdata, n_data, STATUS)
           if (status.ne.adam__ok) then
              write(6,*) 'Error: ',status,' reading Error values from file'
	      call compress(keepscan,cline,c)
              go to 998
           elseif (mod(n_data,2).ne.0) then
              write(6,*) 'Missing Data - ON/OFF mismatch'
              write(6,*) 'Trying to continue'
           endif

*           write(6,*) 'Error: ',status,' Values:'
*           write(6,'(3(e20.10,4x))') (errdata(i),i=1,n_data)
*           write(6,*)
*           write(20,*) 'Error: ',status,' Values:'
*           write(20,'(3(e20.10,4x))') (errdata(i),i=1,n_data)
*           write(20,*)

*   April 27, 1992 (JACH::FJO)
*
*   REMOVE SPIKED PAIR(S)
*
*   Check errdata pairs for SPIKED data, i.e. errdata(i or i+1 ) < -1.0 E+38.
*   Remove spiked pair. (For no cycle reversal, remaining data pairs are
*   still in phase; but for cycle reversal all remaining data pairs must
*   be re-phased by changing signs.)
*   For spiked pair in the data, only the remaining (downstream) data pairs
*   are affected.
* 
*   This method will also remove unrecorded data from an aborted scan
*   since the errdata(i) field is identical for a spike as for an aborted
*   scan.
*
*   SIGNAL VALUE is On-Off by definition.
*
*   NON CYCLE REVERSAL data:
*
*   Spiked--+
*           |
*           V
*   On,Off  On,Off  On,Off  On,Off  On,Off
*   A1,B1,  A2,B2,  A3,B3,  A4,B4,  A5,B5
*       
*
*                 Spiked--+
*                         |
*                         V
*   SIGNAL VALUES: A1-B1, A2-B2, A3-B3, A4-B4, A5-B5
*   
*   SIGNAL VALUES (after removal of spike data pair):
*
*                  A1-B1, A3-B3, A4-B4, A5-B5
*
*   CONCLUSION: for non cycle reversal, when a data pair is removed
*               the remaining data pairs are sign wise still in phase.
*
*   CYCLE REVERSAL data:
*
*  Spiked---+---------------+
*           |               |
*           V               V
*   On,Off  Off,On  On,Off  Off,On  On,Off
*   A1,B1,  A2,B2,  A3,B3,  A4,B4,  A5,B5
*       
*
*   SIGNAL VALUES with out spike data removal:
*  Spiked---+--------------+
*           |              |
*           V              V
*   A1-B1, -A2+B2, A3-B3, -A4+B4, A5-B5
*   
*   Data after removal of first spiked pair (A2,B2) and sign change
*
*             Spiked---+
*                      |
*                      V
*   On,Off  Off,On    On,Off,   Off,On
*   A1,B1,  -A3,-B3,  -A4,-B4,  -A5,-B5
*
*   Data after removal of 2nd spiked pair (-A4,-B4) and sign change
*
*             Spiked---+
*                      |
*                      V
*   On,Off  Off,On    Off,On
*   A1,B1,  -A3,-B3,  -(-A5),-(-B5)


*   SIGNAL VALUES (after removal of all spike data pairs):
*
*                  +On  - Off, -Off  +On,    +On      -Off
*                  +(A1)-(B1), -(-A3)+(-B3), +(-(-A5))-(-(-B5))
*
*     or           A1-B1, A3-B3, A5-B5
*
*   CONCLUSION: for non cycle reversal, when a data pair is removed
*               the remaining data pairs are sign wise still in phase.

           i = 1
           j = 1
           psign = +1
           
           DO WHILE (i.lt.n_data)

              IF ((errdata(i).lt.-1.0e+38).or.(errdata(i+1).lt.-1.0e+38)) then
                 write(6,*) 'Removing Spike datum pair ', (i-1)/2 + 1
                 if (cycle_reversal) psign = -psign
              ELSE
                 errdata(j) = errdata(i)
                 errdata(j+1) = errdata(i+1)
                 IF (psign.gt.0) then
                    data(j) = data(i)
                    data(j+1) = data(i+1)
                 ELSE
                    data(j) = -data(i)
                    data(j+1) = -data(i+1)
                 ENDIF
                 j = j + 2
              ENDIF
              i = i + 2
           ENDDO


* Check for excessive data
           IF ((end+j-1).gt.3000) then
              write(6,*) 'Buffer full.  Rejecting current input file.'
              goto 999
           ENDIF

* Update n_data in case there was a spike
           n_data = j - 1
           nitems = n_data/2

           IF (ans(1:1) .eq. 'y' .or. ans(1:1) .eq. 'Y' .or.
     &         c .le. cline) then
              start = end + 1
              end = end + 2*nitems
           ELSE
              start = 1
              end = 2*nitems
           END if
                      

        indivobs(c)=nitems

        overheads=3.0*(2.0*nitems)      ! 3 secs overhead per nod

        lststart=lst*15.0		! convert hours to degrees
        lstend=lst+(((real(totalint) + overheads)/3600.0)*0.997269) 
      	lstend=lstend*15.0

! 0.997269 lst to solar days

      z1 = 90.0 - el
      h = lstend - ra
      trig = (sind(dec)*sind(lat)) + (cosd(dec)*cosd(h)*cosd(lat))
      z2 = 90.0 - asind(trig)

      amstart(c) = 1.0/cosd(z1)
      amend(c) = 1.0/cosd(z2)


* The normal observing mode is CYCLE_REVERSAL = TRUE
* Successive statements assume this.
* Therefore, for the case CYCLE_REVERSAL = FALSE,
* if the previous data set ended in an odd pair we convert the new data
* by negating odd new data pairs 1, 3, 5, ....
* But if the previous data set ended with and even pair we convert the new data
* by negating even new data pairs 2, 4, 6, ....
*
* In summary, the algorithm to join a new data set to a previous data set
* (where previous data set is assumed CYCLE_REVERSAL=TRUE) is:
*
*       1. If new data set is CYCLE_REVERSAL = TRUE
*          a. If previous data set ends with an even pair
*             then no change to new data points before joining
*          b. If previous data set ends with an odd pair
*             then negate all new data points before joining
*       2. If new data set is CYCLE_REVERSAL = FALSE
*          a. If previous data set ends with an even pair
*             then negate even new data points (2,4,6,...) before joining
*          b. If previous data set ends with an odd pair
*             then negate odd new data points (1,3,5,...) before joining.
*
*  This algorithm is implemented in a more efficient manner as follows:
*
        end_even = mod((start-1)/2,2) .eq. 0
        k=0
        DO i=start, end
           k=k+1
           IF (cycle_reversal) then

              if (end_even) then
                 a(i) = data(k)
              else
                 a(i) = -data(k)
              endif

           ELSE  !(.not.cycle_reversal)

*  Note: new_even is whether the current pair from the new data is even or odd

              new_even = mod((k-1)/2+1,2) .eq. 0

              IF ((end_even .and. new_even) .or.
     :            (.not.(end_even .or. new_even))) then
                 a(i) = -data(k)
              ELSE
                 a(i) = data(k)
              ENDIF

           ENDIF
           b(i) = errdata(k)
        ENDDO

998     IF (c .lt. cline) then       ! more in buffer
          goto 995
        ENDIF
        write(6,*) ' '
        write(6,*) ' '
        write(6,'(A)') '$ coadd more data (Y/N) ? : [N] '
        read(5,'(A)') ans
        IF (ans(1:1) .eq. 'y' .or. ans(1:1) .eq. 'Y') goto 990                

        IF (c .eq. 0) then
          write(6,*) ' No scans, exiting...'
          goto 199
        ENDIF

999     nocoadds = c
        write(6,*) ' '
      	write(6,*) ' total number of coadds : ', c
        counter = end
        nobs = counter/2

	savenobs = nobs !        save nobs for full repeated reduction 

	N=0             ! INITIAL VALUES FOR COUNTERS
	M=0             !

        DO I=1,COUNTER,2
          N=N+1
          MVOLTS(N) = A(I) - A(I+1)                   ! O/P = L - R
        END DO

         write(6,*) ' '
      	 write(6,*) ' total number of paired observations : ', nobs
         write(6,*) ' '
         write(50,*) ' '
      	 write(50,*) ' total number of paired observations : ', nobs
         write(50,*) ' '


        DO I= 1,COUNTER,2
            M=M+1
            ERR(M) = SQRT(B(I)**2 + B(I+1)**2)         ! QUADRATURE
        END DO

        FIRSTERROR = ERR(1)  ! TO BE RETURNED TO SUBROUTINE STDERR

        OPEN (UNIT=11,FILE='mv.dat',STATUS='unknown')

        DO  I=1,NOBS,2
           WRITE (11,*) I, MVOLTS(I), ERR(I)           ! TEMP. FILE FOR011.DAT
           WRITE (11,*) I+1, -1*MVOLTS(I+1), ERR(I+1)  ! ACCOUNTS FOR THE SIGN 
                                                  ! OF EVERY SECOND L-R PAIR
                                                  ! BEING REVERSED
        END DO

	CLOSE (UNIT=11)

101     datafile=11
	CALL MEAN(MVOLTS_ZERO,ERR_ZERO,MEANFIN,MEANMV,ERRORFIN,
     :            MEANERROR,NOBS,OBS,DATAFILE,option)


        call median(mvolts_zero,nobs,mv_median,option)

        CALL STDERR(MEANFIN,MVOLTS,NOBS,MEANMV,FIRSTERROR,STDEV,
     :              ERRORSTDEV,ERRORFIN,NOBS,
     :              OBS,DATAFILE,MV_MEDIAN,option, stderror, signoise)

        CALL TAUCORRECTION(NOBS,MVOLTS,ERR,DATAFILE,indivobs,
     :              amstart,amend,nocoadds,keepscan,filter,option)

	CALL MEAN(MVOLTS_ZERO,ERR_ZERO,MEANFIN,MEANMV,ERRORFIN,
     :              MEANERROR,NOBS,OBS,DATAFILE,option)

* saves the meanvalue of the sample after extinction correction  to allow the
* correct mean to be plotted if the user chooses to repeat the despiking

	taumean = meanfin
*

        call median(mvolts_zero,nobs,mv_median,option)

        CALL STDERR(MEANFIN,MVOLTS,NOBS,MEANMV,FIRSTERROR,STDEV,
     :              ERRORSTDEV,ERRORFIN,NOBS,
     :              OBS,DATAFILE,MV_MEDIAN,option, stderror, signoise)

* similar for sample standard deviation 

	taustdev = stdev
*
	write(6,'(A)') '$ Display data and despike (Y/N) ? : [Y] '
	read(5,'(A)') dans
	if (dans(1:1) .eq. 'N' .or. dans(1:1) .eq. 'n') then 
	  ksnobs = savenobs 
	  goto 103
        end if

102        CALL PLOTGRAPH(NOBS,MVOLTS_ZERO,MVOLTS,ERR_ZERO,ERR,
     :                  DELPAIR,NODEL,MEANFIN,STDEV,DATAFILE,option)

           CALL MEAN(MVOLTS_ZERO,ERR_ZERO,MEANFIN,MEANMV,ERRORFIN,
     :              MEANERROR,NOBS,OBS,DATAFILE,option)

           call median(mvolts_zero,nobs,mv_median,option)

           CALL STDERR(MEANFIN,MVOLTS,NOBS,MEANMV,FIRSTERROR,STDEV,
     :              ERRORSTDEV,ERRORFIN,NOBS,
     :              OBS,DATAFILE,MV_MEDIAN,option, stderror, signoise)

* save the number of samples remaining after despiking so the KS-test
* can be repeated

            ksnobs = nobs
           
c       at this point the data are in the form (pair, mv, err), have been
c       corrected for extinction and have been despiked to a satisfactory 
c       level. now the data are to be divided into a number of 
c       subsamples, KS-tested and coadded. 

        
 103    write(6,'(A)') '$ DO A K-S TEST (Y/N) ? : [Y] '
        read(5,'(A)') kans
        if (kans(1:1) .eq. 'N' .or. kans(1:1) .eq. 'n') goto 110

 104    CALL KSTEST(OBS,MVOLTS_ZERO,ERR_ZERO,DATAFILE,NOBS,option)

        CALL MEAN(MVOLTS_ZERO,ERR_ZERO,MEANFIN,MEANMV,ERRORFIN,
     :              MEANERROR,NOBS,OBS,DATAFILE,option)

        call median(mvolts_zero,nobs,mv_median,option)


        CALL STDERR(MEANFIN,MVOLTS,NOBS,MEANMV,FIRSTERROR,STDEV,
     :              ERRORSTDEV,ERRORFIN,NOBS,
     :              OBS,DATAFILE,MV_MEDIAN,option, stderror, signoise)


 110    optionchar = ' '
        write(6,*) ' '
        write(6,*) ' '        
        write(6,*)' Repeat full reduction with original raw data      : [ 1 ] '
        write(6,*) ' '
        write(6,*)' Repeat despiking of data                          : [ 2 ] '
        write(6,*) ' '
        write(6,*)' Repeat K-S test and selection of subsamples       : [ 3 ] '
        write(6,*) ' '
        write(6,*)' Add new scan(s) and repeat full reduction         : [ 4 ] '
        write(6,*) ' '
        write(6,'(A)') 
     &           '$                                   OPTION [CONTINUE] : '
        read(5,'(a)') optionchar
        if (optionchar(1:1) .eq. ' ') then
          option = 0
        else
          read(optionchar,*) option
        end if

        if (option .ge. 1 .and. option .le. '4') then

          if (option .eq. 2) then
            datafile=14
            write(6,*) ' '
            write(6,*) '... repeat despiking of data ...'
            nobs = savenobs
            meanfin = taumean
            stdev = taustdev
            write(6,*) ' '        
            goto 102
         
          else if (option .eq. 3) then
            datafile=16
            nobs=ksnobs
            write(6,*) ' '
            write(6,*) '... repeating coadding and K-S testing of subsamples ...'
            write(6,*) ' '        
            goto 104
         
          else if (option .eq. 4) then
            datafile=11
            write(6,*) ' '
            write(6,*) '  ADD NEW SCAN(S):'
            nobs = savenobs
            goto 990
         
          else
            datafile=11
            write(6,*) ' '
            write(6,*) '... repeating reduction of original raw data ...'
            nobs = savenobs
            write(6,*) ' '        
            goto 101
         
          end if
        end if

        write(6,*) '  '
        write(6,'(A)') '    Calibrate final reduced signal level (mV) via '
        write(6,'(A)') '$        known responsivity (Jy/mV)     (Y/N) ? : [N] '
        read(5,'(A)') ans

        if (ans(1:1) .eq. 'Y' .or. ans(1:1) .eq. 'y') then
          write(6,*) ' '
          write(6,'(a)') '$ Responsivity : '
          read(5,*) respons

          calmean = meanfin * respons
          calerror = stderror * respons


9040    format(6x,A9,7x,A11,5x,A3)
9050    format(2x,f12.6,1x,a3,2x,f12.6,2x,f7.3)
        WRITE(6,*) '      '
        WRITE(6,9040) 'MEAN (Jy)','err (Jy)','S/N'
        WRITE(6,9050)  calmean,'+/-',calerror,signoise
        WRITE(50,9040) 'MEAN (Jy)','err (Jy)','S/N'
        WRITE(50,9050)  calmean,'+/-',calerror,signoise
        WRITE(6,*) '      '

        end if          

!--------------------------------------------------------------------------!
!                                CLEAR UP                                  !
!--------------------------------------------------------------------------!

        write(6,*) '  '
        write(6,*) '  '
        write(6,'(a26,2x,a20,a4)') ' Summary of reduction in ',summary,'.DAT'
        write(6,*) '  '

 199    write(6,*) '  '
        write(6,*) ' ... tidying up . '

        close(unit=11,status='delete')
        close(unit=12,status='delete')
        close(unit=13,status='delete')
        close(unit=14,status='delete')
        close(unit=16,status='delete')
        if (datafile .eq. 20) then
          close(unit=20,status='delete')
        endif

        END

!--------------------------------------------------------------------------!
!                       SUBROUTINE PARSELINE                               !
!--------------------------------------------------------------------------!

       SUBROUTINE PARSELINE(keepscan,cline)

c  Get scans from user or file. Incorporates simple line parser.
c  Not elegant, but works.

       IMPLICIT NONE

      integer keepscan(100)
      integer cline

      integer      i, j, k
      integer      iscan, fscan, lscan
      integer      range
      character*4  scan
      character*40 filename
      character*81 line
      logical      infile

      integer*4    bi
      character*4  bc
      equivalence(bi,bc)

      bi = 0
      infile = .false.

      write(6,'(A)') 
     & '$  Give UKT14 GSD scan number(s) (e.g. 1-10,12 14,15 20:25) or'
      write(6,'(A)') 
     & '$  Name of file with numbers (same format, multiple lines allowed).'

 5001 continue
      line(1:80) = ' '
      if (.not. infile) then
        write(6,'(A)') 
     &     '$  #(s) or Name [CONTINUE]: '
        read(5,'(A80)') line
      else   
         read(9,'(A80)',end=5090) line
      endif

      i = 1
      do while (line(i:i) .eq. ' ' .and. i .lt. 81)
        i = i + 1
      enddo
      j = 80
      do while (line(j:j) .eq. ' ' .and. j .ge. 1)
        j = j - 1
      enddo

      if (i .eq. 81) goto 5099                 ! Blank line: finished

c (VAX:) bc(1:1) = line(i:i)
      bc(4:4) = line(i:i)
      if ( bi .lt. 48 .or. bi .gt. 57) then    ! Starts with char: filename
        filename(1:40) = ' '
        if ( (j-i) .gt. 40 ) then
          write(6,*)' '
          write(6,'(A)') 
     & '    ERROR: Filename too long (max. 40 characters)'
          write(6,*)' '
          goto 5001
        endif

        infile = .true.
        filename = line(i:j)
        open(unit=9, file=filename, status='old')
        line(1:80) = ' '
        read(9,'(A80)',end=5090) line
        do while (line(i:i) .eq. ' ' .and. i .lt. 81)
          i = i + 1
        enddo
        do while (line(j:j) .eq. ' ' .and. j .ge. 1)
          j = j - 1
        enddo
        if (i .eq.81) goto 5099                  ! Blank line: finished
      endif

      j = j + 1                                  ! Add terminating blank

      k = i
      range = 0
      do while (k .le. j)

        if (line(k:k) .ne. ',' .and. line(k:k) .ne. '-' 
     &      .and. line(k:k) .ne. ':' .and. line(k:k) .ne. ' ') then
          k = k + 1
        else
          if (range .eq. 0) then
            read(line(i:(k-1)),'(I)') fscan
            lscan = fscan
          else     
            read(line(i:(k-1)),'(I)') lscan
          endif

          if (line(k:k) .eq. ',' .or. line(k:k) .eq. ' ') then
            do iscan = fscan, lscan
               cline = cline + 1
               if (cline .gt.100) goto 5095
               keepscan(cline) = iscan 
            enddo
            range = 0
          else if (line(k:k) .eq. '-' .or. line(k:k) .eq. ':') then
            range = 1
          endif

          k = k + 1
          i = k

        endif

      enddo
c
c Return for next line
c
      goto 5001

 5090 continue
      close(unit=9)
      infile = .false.
      goto 5001

 5095 continue
      write(6,*)' '
      write(6,'(''    ERROR: max scans (100) exceeded at scan '',I4)')
      write(6,*)' '

 5099 continue  
      if (infile) then
        close(unit=9)
      endif
      return
      end

!--------------------------------------------------------------------------!
!                       SUBROUTINE COMPRESS                                !
!--------------------------------------------------------------------------!

       SUBROUTINE COMPRESS(keepscan,cline,c)

c  Eliminate faulty scan from scan list

       IMPLICIT NONE

      integer keepscan(100)
      integer cline,c,i

        write(6,'(''    Ignoring scan '',I4,''; Eliminated from list.'')')
     &               keepscan(c)
      do i = c, cline-1
        keepscan(i) = keepscan(i+1)
      enddo
      cline = cline - 1
      c = c - 1

      return
      end

!--------------------------------------------------------------------------!
!                       SUBROUTINE TAUCORRECTION                           !
!--------------------------------------------------------------------------!

       SUBROUTINE TAUCORRECTION(NOBS,MVOLTS,ERR,DATAFILE,indivobs,
     :             amstart,amend,nocoadds,keepscan,filter,option)

c  correction for atmospheric attenuation. 
c  single or multiple opacities allowed for coadded data.
c  linear interpolation between start and finish airmass.
 

                         IMPLICIT NONE

        REAL MVOLTS_ZERO(3000),MVOLTS(3000),ERR(3000),TAU,AIRMASS
        REAL AMSTART(100),AMEND(100),AMINC, OBS(3000), ERR_ZERO(3000)
        real tausingle,taustart(100),tauend(100),tauinc     
        INTEGER NOBS,I,DATAFILE,nocoadds,indivobs(100),j,k,option

        integer keepscan(100)
        CHARACTER*4 ans,answer
        character*8 filter

3001   format(1x,a6,1x,i4,2x,a17,1x,f10.4,1x,a3,1x,f10.4)
3002   format(1x,i4,3x,f10.4,5x,f6.3,4x,f11.5,3x,f11.5)
3003   format(3x,a1,8x,a7,6x,a3,11x,a2,5x,a19)
3004   format(3x,a11,1x,a4,1x,a4)
3005   format(3x,a11,1x,a4,1x,a30)

        WRITE(6,'(A)') '$ Correction for atmospheric extinction (Y/N) ? : [N] '
        read(5,'(A)') answer
        if (answer(1:1) .eq. 'Y' .or. answer(1:1) .eq. 'y') then
          write(6,*) ' '
          write(6,'(A)') '$ Do you wish a single value of tau (Y/N) ? : [Y] '
            read(5,'(A)') ans
            if (ans(1:1) .eq. 'N' .or. ans(1:1) .eq. 'n') then
              write(6,*) ' '
              write(6,3005) 'opacity at ',filter,'mm. for airmasses 
     : given below '               
              write(50,3005) ' opacity at ',filter,'mm. for airmasses 
     : given below '               
              write(6,*) ' '
              do j=1,nocoadds
                write(6,*) 'airmass ',amstart(j), amend(j)
                read(5,*)  taustart(j), tauend(j)
                write(6,*) ' '
                write(50,*) 'airmass ',amstart(j), amend(j)
                write(50,*) '  tau   ',taustart(j), tauend(j)
               end do
            else
              write(6,*) ' '
              write(6,3004) ' opacity at ',filter,'mm.'    
              read(5,*) tausingle
              write(50,3004) ' opacity at ',filter,'mm. '    
              write(50,*) tausingle
            end if
        end if

       OPEN(UNIT=11,FILE='mv.dat',STATUS='OLD')
       OPEN(UNIT=14,FILE='mv_taucorrected.dat',STATUS='unknown')

* reset counter for scan number 
       i=0

       do j=1,nocoadds
       AMINC = (AMEND(j)-AMSTART(j))/REAL(indivobs(j))
       tauinc = (tauend(j)-taustart(j))/real(indivobs(j))
       print*, ' '
       print*, ' '
       print*, ' '
       write(6,3001) ' scan ',keepscan(j), ' airmass range = ', 
     :                amstart(j),' - ', amend(j)
       print*, ' '
       write(6,3003) '#', 'airmass', 'tau', 'mV', 'mV (tau corrected)'
       print*, ' '

       DO k=1,indivobs(j)
         i=i+1
         READ(11,*) OBS(I), MVOLTS(I), ERR(I)
         AIRMASS=AMSTART(j)+((k-1)*AMINC)

           if (answer(1:1) .eq. 'Y' .or. answer(1:1) .eq. 'y') then
            if (ans(1:1) .eq. 'N' .or. ans(1:1) .eq. 'n') then
             tau = taustart(j)+((k-1)*tauinc)
           else
             tau = tausingle
           end if
         else
* no extinction correction
           tau = 0.0             
         end if

         MVOLTS_ZERO(I)=MVOLTS(I)*EXP(TAU*AIRMASS)         
         ERR_ZERO(I)=ERR(I)*EXP(TAU*AIRMASS)
         WRITE(14,*) OBS(I), MVOLTS_ZERO(I), ERR_ZERO(I)
      
          write(6,3002) i, airmass, tau, mvolts(i), mvolts_zero(i)

       END DO
       end do

       print*,' '
       print*,' '
       print*,' Data corrected for atmospheric attenuation'
       print*,' '

        datafile=14    ! allows subroutines mean,median,stderr to
                       ! recalulate         
 
       CLOSE(UNIT=11)
       CLOSE(UNIT=14)
       END

!-----------------------------------------------------------------------!
!                      SUBROUTINE PLOTGRAPH                             !
!-----------------------------------------------------------------------!
 
        SUBROUTINE PLOTGRAPH(NOBS,MVOLTS_ZERO,MVOLTS,ERR_ZERO,ERR,
     :                   DELPAIR,NODEL,MEANFIN,STDEV,DATAFILE,option)
c       Graphics by PGPLOT
c       plots signal and error against pair no. repeats it with tau 
c       corrected signals and also automatic or interactive despiking.

                  IMPLICIT NONE

        REAL MEANFIN,STDEV,NEWOBS,ERRORFIN,ERRORSTDEV
        REAL PAIR(3000), MVMIN, MVMAX, ERRMIN, ERRMAX
        REAL YMIN, YMAX, MVOLTS_ZERO(3000), MVOLTS(3000)
        REAL ERR_ZERO(3000), ERR(3000), X,Y,DESPAIR(3000) 
        REAL MVOLTS_DESP(3000),ERR_DESP(3000),SIG3,PLUS,MINUS

	real piksrt(4), pik

        INTEGER I,II,J,JJ,K,KK,NOBS,DELPAIR(300),NODEL,DATAFILE
	INTEGER SAMPLESIZE, OLDNOBS, n, m,option

        CHARACTER*1 CH
        CHARACTER*5 ANS, DANS
        CHARACTER*20 DEVICE, HDEVICE, ODEVICE
        LOGICAL     HCOPY

        OPEN(UNIT=11,FILE='mv.dat',STATUS='OLD')
        OPEN(UNIT=14,FILE='mv_taucorrected.dat',STATUS='OLD')

        MVMIN=1e37
        MVMAX=-1e37
        ERRMIN=ERR_ZERO(1)
        ERRMAX=ERR_ZERO(1)
        HCOPY = .FALSE.
        HDEVICE = '/PS'

        DO I=1,NOBS
            PAIR(I)=REAL(I)
            READ(11,*) JJ, MVOLTS(I), ERR(I)
            READ(14,*) II, MVOLTS_ZERO(I), ERR_ZERO(I)

            IF (MVOLTS_ZERO(I) .LT. MVMIN) THEN
                MVMIN=MVOLTS_ZERO(I)
            END IF     
            IF (MVOLTS_ZERO(I) .GT. MVMAX) THEN
                MVMAX=MVOLTS_ZERO(I)  
            END IF     
        END DO

        SIG3 = 3.0 * STDEV
        PLUS = MEANFIN + SIG3
        MINUS = MEANFIN - SIG3

	piksrt(1) = abs(mvmax)   ! take absolute value of 3 sig. limits
	piksrt(2) = abs(mvmin)   ! and upper/lower data and sort to
	piksrt(3) = abs(plus)    ! find max deviation from zero.
	piksrt(4) = abs(minus)   !

        call sort(4,piksrt)      !

	ymax =  1.15*piksrt(4)    ! max, min values for plotting range
	ymin = -1.15*piksrt(4)    !

        CLOSE(UNIT=11)
        CLOSE(UNIT=14)

 666    CONTINUE
        write(6,*) '   '
        WRITE(6,'(A50,$)') ' Enter plot device [/TEK,/XSERVE,/PS,etc] '
        READ (5,'(A)') DEVICE

*       List all PGPLOT devices
        IF (DEVICE(1:1) .EQ. '?') THEN
           CALL PGLDEV
           GOTO 666
        ENDIF

        ODEVICE = DEVICE
 700    CALL UUCASE(DEVICE)
        PRINT *, 'Device is ',DEVICE
        IF (DEVICE(2:4).EQ.'TEK') THEN
           WRITE(6,*) 'Please note, to continue from a plot'
           WRITE(6,*) 'Simply press the Return key'
           CALL SLEEP(2)             !PAUSE           
           CALL MODE401X
        ENDIF
        CALL PGBEGIN(0,DEVICE,2,1)
        CALL PGENV(0.,REAL(NOBS)+1.0,YMIN,YMAX,0,0)
        CALL PGLABEL('pair no. (L-R)', 'Signal - measured 
     :               (millivolts)','JCMT/UKT14 data')
        CALL PGPOINT(NOBS,PAIR,MVOLTS,6)
        CALL PGMOVE(0.,0.)
        CALL PGDRAW(REAL(NOBS)+1.0,0.)

        CALL PGENV(0.,REAL(NOBS)+1.0,YMIN,YMAX,0,0)
        CALL PGLABEL('pair no. (L-R)', 'Signal - corrected 
     :               (millivolts)','JCMT/UKT14 data')
        CALL PGPOINT(NOBS,PAIR,MVOLTS_ZERO,6)
        CALL PGMOVE(0.,0.)
        CALL PGDRAW(REAL(NOBS)+1.0,0.)
        CALL PGSLS(2)
        CALL PGMOVE(0.,MEANFIN)
        CALL PGDRAW(REAL(NOBS)+1.0,MEANFIN)
        CALL PGSLS(4)
        CALL PGMOVE(0.,PLUS)
        CALL PGDRAW(REAL(NOBS)+1.0,PLUS)
        CALL PGMOVE(0.,MINUS)
        CALL PGDRAW(REAL(NOBS)+1.0,MINUS)

!-------------------------
c temporary patch for counting problem, EEE reset in AUTODESPIKE

        oldnobs=nobs
!-------------------------

        IF (DEVICE(1:3).EQ.'TEK') THEN
           READ (5,'(A)') ANS
           CALL MODEVT300
        ENDIF
        print*,' '
	print*,' '
	print*,' MEAN  shown as a dashed line'
        print*,' 3 SIGMA LIMITS shown as dotted lines'
	print*,' '

        IF (HCOPY) THEN 
            CALL PGEND
            goto 950
        ENDIF

        WRITE(6,*) ' '
        WRITE(6,'(A)') '$ DO YOU WANT TO DESPIKE (Y/N) ? : [Y] '

        READ(5,'(A)') DANS
        IF (DANS(1:1) .EQ. 'N' .OR. DANS(1:1) .EQ. 'n') THEN
            DATAFILE=14				! tau-corrected file
            CALL PGEND
            goto 910
        ELSE
	    DATAFILE=16				! despiked data file
            WRITE(6,'(A)') '$ MANUAL OR AUTO (M/A) - [A]: '
            READ(5,'(A)') ANS
            IF (ANS(1:1) .EQ. 'M' .OR. ANS(1:1) .EQ. 'm') THEN

               IF (DEVICE(1:3).EQ.'TEK') CALL MODE401X
               CALL PGEND
               CALL PGBEGIN(0,DEVICE,1,1)

               CALL PGENV(0.,REAL(NOBS)+1.0,YMIN,YMAX,0,0)
               CALL PGLABEL('pair no. (L-R)', 'Signal - corrected 
     :               (millivolts)','JCMT/UKT14 data')
               CALL PGPOINT(NOBS,PAIR,MVOLTS_ZERO,6)
               CALL PGMOVE(0.,0.)
               CALL PGDRAW(REAL(NOBS)+1.0,0.)

                   J=0
800                CALL PGCURSE(X,Y,CH)
                   J=J+1
                   DELPAIR(J)=JNINT(X)    

               write(6,*) '   '
               WRITE(6,*) ' data at coords (pair, mV) : ',DELPAIR(J),Y
               write(6,*) ' '
               WRITE(6,'(A)') '$  ... again? '
               READ(5,'(A)') ANS
               IF (ANS(1:1) .EQ. 'Y' .OR. ANS(1:1) .EQ. 'y') GOTO 800
               NODEL = J      ! no. of data deleted
               write(6,*) '  '
               write(6,*) ' pairs of data deleted : ',nodel
               CALL MANDESPIKE(MVOLTS_ZERO,ERR_ZERO,NOBS,DELPAIR,
     :                       NODEL,SAMPLESIZE,option)       
            ELSE
               IF (DEVICE(1:3).EQ.'TEK') CALL MODE401X
               CALL AUTODESPIKE(MEANFIN,STDEV,MVOLTS_ZERO,ERR_ZERO,
     :                       NOBS,oldnobs,SAMPLESIZE,option)
            END IF

            CALL PGEND
            IF (DEVICE(1:3).EQ.'TEK') THEN
               READ (5,'(A)') ANS
               CALL MODEVT300
            ENDIF
        END IF


	write(6,*) ' '
        WRITE(6,'(A)') '$ DO YOU WANT TO RE-PLOT (Y/N) ? : [N] '
        READ(5,'(A)') ANS
        IF (ANS(1:1) .EQ. 'Y' .OR. ANS(1:1) .EQ. 'y') THEN
           OPEN(UNIT=16,FILE='despiked.dat',STATUS='OLD')
           K=0
           DO I=1,NOBS
           READ(16,*,END=900) KK,MVOLTS_DESP(I),ERR_DESP(I)
           DESPAIR(I)=KK
           K=K+1
           END DO
   
900        write(6,*) '   '
           IF (DEVICE(1:3).EQ.'TEK') CALL MODE401X

           CALL PGBEGIN(0,DEVICE,2,1)
           CALL PGENV(0.,REAL(OLDNOBS)+1.0,YMIN,YMAX,0,0)
           CALL PGLABEL('pair no. (L-R)', 'Signal - tau corrected
     :               (millivolts)','JCMT/UKT14 data')
           CALL PGPOINT(OLDNOBS,PAIR,MVOLTS_ZERO,6)
           CALL PGMOVE(0.,0.)
           CALL PGDRAW(REAL(OLDNOBS)+1.0,0.)

           CALL PGENV(0.,REAL(OLDNOBS)+1.0,YMIN,YMAX,0,0)
           CALL PGLABEL('pair no. (L-R)', 'Signal - tau corrected 
     :               (millivolts)','JCMT/UKT14 despiked data')
           CALL PGPOINT(K,DESPAIR,MVOLTS_DESP,6)
           CALL PGMOVE(0.,0.)
           CALL PGDRAW(REAL(OLDNOBS)+1.0,0.)
           CALL PGEND
           IF (HCOPY) GOTO 940
           IF (DEVICE(1:3).EQ.'TEK') THEN
               READ (5,'(A)') ANS
               CALL MODEVT300
           ENDIF
        END IF

       print*,' '
       print*,' '
       print*,' Data corrected for atmospheric attenuation and despiked'
       print*,' '

 910   write(6,*) ' '
        WRITE(6,'(A)') '$ DO YOU WANT A HARDCOPY (Y/N) ? : [N] '
        READ(5,'(A)') ANS     
        IF (ANS(1:1) .EQ. 'Y' .OR. ANS(1:1) .EQ. 'y') THEN
           HCOPY = .TRUE.
           DEVICE = HDEVICE
           IF (DANS(1:1) .EQ. 'N' .OR. DANS(1:1) .EQ. 'n') THEN
             GOTO 700
           ELSE
             GOTO 900
           ENDIF
        ENDIF

 940    CLOSE(UNIT=16)  
 950    DEVICE = ODEVICE
        END

* Routines to switch from Tektronix to VT300 terminal modes.
* Terminal must be vt330

      subroutine mode401x
*  enter VT330 tektronix mode.
      character esc
      byte esc1
      equivalence (esc1,esc)
      data esc1 /27/
      write(*,'(a1,a1,a5)') ' ', esc, '[?38h'
      end
      subroutine modevt300
*  enter VT330 tektronix mode.
      character esc
      byte esc1
      equivalence (esc1,esc)
      data esc1 /27/
*  that is a lower case L at the end of [?38l
      write(*,'(a1,a1,a5)') ' ', esc, '[?38l'
      end
*

*--------------------------------------------------------------------------!
*                            SUBROUTINE  MEAN                              !
*--------------------------------------------------------------------------!

       	SUBROUTINE MEAN(YY,ZZ,MEANFIN,MEANMV,ERRORFIN,MEANERROR,
     :                  NOBS,OBS,DATAFILE,option)

*           CALCULATES MEAN DATA AND MEAN ERROR FROM INDIVIDUAL L-R's'

* YY           REAL ARRAY(3000) MILLIVOLTS FROM BOLOMETER
* ZZ           REAL ARRAY(3000) ERROR CALCULATED IN QUADRATURE FROM L & R
* MEANFIN      REAL            FINAL MEAN FROM FULL OBSERVATION TO BE RETURNED 
*                               TO SUBROUTINE STDERR
* MEANMV       REAL ARRAY(3000) RUNNING MEAN RETURNED TO SUBROUTINE STDERR
*                               TO CALCULATE RUNNING GAIN IN SIG/NOISE.
* ERRORFIN     REAL             FINAL MEAN ERROR FROM SAMPLE TO BE RETURNED TO
*                               SUBROUTINE STDERR
* MEANERROR    REAL ARRAY(3000) RUNNING MEAN ERROR RETURNED 
*                               TO SUBROUTINE STDERR
* NOBS         INTEGER          NO. OF PAIRS IN OBSERVATION
* OBS          REAL ARRAY(3000) HANGS ONTO IDENTIFICATION OF PAIR

        IMPLICIT NONE

        REAL YY(3000), ZZ(3000), MEANFIN, ERRORFIN
	REAL OBS(3000), MEANMV(3000), ERROR(3000)
        REAL SN(3000), STDEV(3000)
	REAL SUMMVOLTS, WTSUMMVOLTS 
        REAL SUMERROR, MEANERROR(3000) 

	INTEGER N, NOBS, DATAFILE,option

        IF (DATAFILE .EQ. 11) THEN
           OPEN(UNIT=11,FILE='mv.dat',STATUS='OLD')
        ELSE IF (DATAFILE .EQ. 14) THEN 
           OPEN(UNIT=14,FILE='mv_taucorrected.dat',
     :          STATUS='OLD')
        ELSE IF (DATAFILE .EQ. 16) THEN
           OPEN(UNIT=16,FILE='despiked.dat',
     :          STATUS='OLD')
        ELSE
           OPEN(UNIT=20,FILE='coadded.dat',STATUS='OLD')
        END IF

        SUMMVOLTS = 0.0                                ! RESET
        SUMERROR = 0.0
        
        DO  N=1,NOBS
           READ(DATAFILE,*,END=2000) OBS(N), YY(N), ZZ(N)
               
           IF (N.NE.1) THEN

              SUMMVOLTS = SUMMVOLTS + YY(N) ! RUNNING SUM OF MVOLTS
              MEANMV(N) = SUMMVOLTS / N     ! RUNNING MEAN OF MVOLTS 
              SUMERROR  = SUMERROR + ZZ(N)  ! RUNNING SUM OF ERRORS
              MEANERROR(N) = SUMERROR / N   ! RUNNING MEAN OF ERRORS

           ELSE
              MEANMV(1) = YY(1)
              SUMMVOLTS = YY(1)
              MEANERROR(1) = ZZ(1)
              SUMERROR = ZZ(1)
           ENDIF               
           WRITE(12,*) OBS(N), MEANMV(N), MEANERROR(N)
        END DO
           
 2000   MEANFIN = MEANMV(N-1)     ! RETURNS FINAL MEAN TO 'MEAN' 
        ERRORFIN = MEANERROR(N-1) !     SIMILAR

        CLOSE (UNIT=11)
        CLOSE (UNIT=12)
        CLOSE (UNIT=14)
        CLOSE (UNIT=16)         
        CLOSE (UNIT=20) 
        END


*-----------------------------------------------------------------------------!
*                            SUBROUTINE MEDIAN                                !
*-----------------------------------------------------------------------------!

      SUBROUTINE MEDIAN(X,N,MV_MEDIAN,option)

      REAL X(3000), MV_MEDIAN
      
      INTEGER N2,N,option

      CALL SORT(N,X)

      N2=N/2
      IF(2*N2.EQ.N)THEN
         MV_MEDIAN=0.5*(X(N2)+X(N2+1))
      ELSE
         MV_MEDIAN=X(N2+1)
      ENDIF    
   
      RETURN
      END

*----------------------------------------------------------------------------!
*                             SUBROUTINE SORT                                !
*----------------------------------------------------------------------------!

      SUBROUTINE SORT(N,RA,option)

c     sorts data for KS-test
      
      REAL RA(3000),RRA

      INTEGER L,N,I,J,IR,option

      L=N/2+1
      IR=N
 10   CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         RRA=RA(L)
      ELSE
         RRA=RA(IR)
         RA(IR)=RA(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
 20   IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
         ENDIF
         IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         ENDIF
         GOTO 20
      ENDIF
      RA(I)=RRA
      GOTO 10
      END


*-----------------------------------------------------------------------------!
*                            SUBROUTINE STDERR                                !
*-----------------------------------------------------------------------------!

      SUBROUTINE STDERR(XX,YY,ZZ,DD,EE,FF,ERRORSTDEV,GG,RR,
     +     OBS,DATAFILE,MV_MEDIAN,option,stderror,signoise)

*                CALCULATES STANDARD ERROR  IN THE MEAN

* XX          REAL              FINAL SAMPLE MEAN RETURNED FROM SUBROUTINE MEAN
* YY          REAL ARRAY(3000)  INDIVIDUAL L-R's' FROM MV.DAT (UNIT 11)
* ZZ          INTEGER           NO. OF OBSERVATIONS FROM TOP LEVEL (NOBS)
* DD          REAL ARRAY(3000)  RUNNING MEANS FROM SUBROUTINE MEAN
* EE          REAL              STDEV OF 1st OBSERVATION, BEHAVES AS STDERROR
*                                    IN 1st CALCULATION OF SIGNAL TO NOISE
* FF          REAL              STANDARD DEVIATION OF SAMPLE
* ERRORSTDEV  REAL              STANDARD DEVIATION OF ERRORS IN SAMPLE
* GG          REAL              FINAL SAMPLE MEAN ERROR FROM SUBROUTINE MEAN
* RR          INTEGER           NO OF OBSERVATIONS
* OBS         REAL ARRAY(3000)  HANGING ONTO IDENTIFICATION NO. OF PAIR

      IMPLICIT NONE

      REAL XX, YY(3000), AA(3000), BB(3000), CC(3000), DD(3000), EE, FF
      REAL SUMSQRESID, ERROR, SIGNOISE, RUNERROR, RUNSTDEV, GG, HH
      REAL RUNSIGNOISE, SQRESID, FIRSTERROR, SUMSQERROR, ERRORSTDEV
      REAL SUMTEMPMEAN, SUMTEMPERR, OBS(3000),stderror
      REAL MV_MEDIAN
      
      INTEGER ZZ, I, COUNTER, M, J, DUMMY, NEWCOUNTER, K, L
      INTEGER DATAFILE, II, RR, option


      IF (DATAFILE .EQ. 11) THEN
         OPEN(UNIT=11,FILE='mv.dat',STATUS='OLD')
      ELSE IF (DATAFILE .EQ. 14) THEN
         OPEN(UNIT=14,FILE='mv_taucorrected.dat',
     :        STATUS='OLD')
      ELSE IF (DATAFILE .EQ. 16) THEN
         OPEN(UNIT=16,FILE='despiked.dat',
     :        STATUS='OLD')
      ELSE 
         OPEN (UNIT=20,FILE='coadded.dat',STATUS='OLD')
      END IF

      OPEN(UNIT=12,FILE='for012.dat',STATUS='unknown')
      OPEN(UNIT=13,FILE='stderr.dat',STATUS='unknown')

      COUNTER = 0               ! RESET COUNTER
      SUMSQRESID = 0.0
      SUMSQERROR = 0.0


      DO  I=1,RR 
         READ(DATAFILE,*,END=3000) AA(I), BB(I), CC(I)
         COUNTER = COUNTER + 1 
         SUMSQRESID = SUMSQRESID + ((BB(I) - XX)**2) ! SUMMED SQUARED RESIDs
                                                     ! OF DATA
         SUMSQERROR = SUMSQERROR + ((CC(I) - GG)**2) ! SUMMED SQUARED RESIDs
                                                        ! OF ERRORS
      END DO

3000    FF  = (SUMSQRESID/(COUNTER-1))**0.5             ! ST.DEV. IN MEAN
        ERRORSTDEV  = (SUMSQERROR/(COUNTER-1))**0.5     ! ST.DEV. IN ERROR
        ERROR = FF/((COUNTER)**0.5)                     ! STANDARD ERROR
	stderror = error

        SIGNOISE = XX / ERROR

3040    format(6x,A9,6x,A6,7x,A8,5x,A3)
3050    format(2x,f12.6,2x,f12.6,2x,f12.6,2x,f7.3)
        WRITE(6,*) '      '
        WRITE(6,3040) 'MEAN (mV)','MEDIAN','err (mV)','S/N'
        WRITE(6,3050)     XX, MV_MEDIAN, ERROR, SIGNOISE
        WRITE(6,*) '      '
        WRITE(6,*) '      '

* output results to 'SUMMARY.DAT'

	if (datafile .eq. 11) then
          WRITE(50,*) '      '
	  write(50,*) '       Raw data - '
	  write(50,*) '       -----------'
          WRITE(50,*) '      '
	else if (datafile .eq. 14) then
          WRITE(50,*) '      '
	  write(50,*) '       extinction-corrected data - '
	  write(50,*) '       ----------------------------'
          WRITE(50,*) '      '
	else if (datafile .eq. 16) then
          WRITE(50,*) '      '
	  write(50,*) '       despiked data - '
	  write(50,*) '       ---------------'
          WRITE(50,*) '      '
	else 
          WRITE(50,*) '      '
	  write(50,*) '       coadded and K-S tested data - '
	  write(50,*) '       ------------------------------'
          WRITE(50,*) '      '
	end if
          WRITE(50,3040) 'MEAN (mV)','MEDIAN','err (mV)','S/N'
          WRITE(50,3050)     XX, MV_MEDIAN, ERROR, SIGNOISE
          WRITE(50,*) '      '
          WRITE(50,*) '      '


!-----------------------------------------------------------------------------!
!        CALCULATES THE RUNNING SIGNAL TO NOISE AND STANDARD ERROR            !
!-----------------------------------------------------------------------------!
        

       DO I=1,3000
          READ(12,*,END=3010) J, DD(I), CC(I)
       END DO


3010   DO I=1, COUNTER
           
         SQRESID=0

         DO J=1,I

           IF (I .EQ. 1 ) THEN
                RUNSIGNOISE = DD(1) / EE              ! FIRST (L-R)/ FIRST ERROR
           ELSE
       	        SQRESID = SQRESID + ((BB(J) - DD(I))**2) ! SUMMED SQUARED RESIDs
                RUNSTDEV = (SQRESID/REAL(I-1))**0.5
                RUNERROR = RUNSTDEV/(REAL(I)**0.5)
                RUNSIGNOISE = DD(I) / RUNERROR
           ENDIF
         END DO
         WRITE(13,*) I, RUNERROR, RUNSIGNOISE
      END DO
        
!-----------------------------------------------------------------------------!
!  CALCULATES THE RMS SEQUENTIALLY THROUGH THE SAMPLE IGNORING ONE DATA PT.   !
!-----------------------------------------------------------------------------!

        NEWCOUNTER = COUNTER

        DO L=1,NEWCOUNTER

          COUNTER = 0                                     ! RESET COUNTER
          DUMMY = 0
          SUMSQRESID = 0.0
          SUMSQERROR = 0.0
          SUMTEMPMEAN = 0.0
          SUMTEMPERR = 0.0
  
          K=K+1    ! ALLOWS THE READ TO IGNORE EVERY ITERATION THE Kth TERM
                   ! WHEN CALCULATING THE RMS OF EACH ITERATION.

!---------------
! CALCULATES THE APPROPRIATE MEANS ( ERROR AND MEAN DATA)

          DO  II=1,NEWCOUNTER
               COUNTER = COUNTER + 1 
             IF (K. NE. COUNTER) THEN
  	       SUMTEMPMEAN = SUMTEMPMEAN + BB(II)
               SUMTEMPERR  = SUMTEMPERR  + CC(II)
             DUMMY = DUMMY + 1     ! NO. OF DATA IN MEAN VALUE     
             ENDIF
          END DO

             XX = SUMTEMPMEAN / DUMMY
             GG = SUMTEMPERR  / DUMMY

             COUNTER = 0
!---------------

          DO  I=1,NEWCOUNTER
               COUNTER = COUNTER + 1 
             IF (K. NE. COUNTER) THEN
  	       SUMSQRESID = SUMSQRESID + ((BB(I) - XX)**2)  ! SUMMED SQUARED 
                                                            ! RESIDs OF DATA
               SUMSQERROR = SUMSQERROR + ((CC(I) - GG)**2)  ! SUMMED SQUARED 
                                                            ! RESIDs OF ERRORS
             ENDIF
          END DO
    
           FF  = (SUMSQRESID/(COUNTER-1))**0.5             ! ST.DEV. IN MEAN
           ERRORSTDEV  = (SUMSQERROR/(COUNTER-1))**0.5     ! ST.DEV. IN ERROR
           ERROR = FF/((COUNTER)**0.5)                     ! STANDARD ERROR
           SIGNOISE = XX / ERROR
        END DO

       CLOSE (UNIT=11)
       CLOSE (UNIT=12)
       CLOSE (UNIT=13)
       CLOSE (UNIT=14)
       CLOSE (UNIT=16)    
       CLOSE (UNIT=20)      
       END


!----------------------------------------------------------------------------!
!                       SUBROUTINE AUTODESPIKE                               !
!----------------------------------------------------------------------------!

!       DESPIKES ALL DATA VALUES GREATER THAN ? SIGMA FROM SAMPLE MEAN

       SUBROUTINE AUTODESPIKE(AAA,BBB,CCC,DDD,EEE,oldnobs,FFF,option)

! AAA   REAL             SAMPLE MEAN
! BBB   REAL             STANDARD DEVIATION OF SAMPLE
! CCC   REAL ARRAY(3000) INDIVIDUAL DATA
! DDD   REAL ARRAY(3000) INDIVIDUAL ERRORS
! EEE   INTEGER          NO OF DATA IN ORIGINAL SAMPLE
! FFF   INTEGER          NO OF DATA IN NEW SAMPLE

	REAL AAA, BBB, CCC(3000), DDD(3000), PLUS, MINUS
        REAL SIGMA

        INTEGER I, EEE, samplesize, FFF, oldnobs, option
	character*4 ans
	character*10 sigchar

        FFF = 0    ! RESET
      
        OPEN (UNIT=16,FILE='despiked.dat',STATUS='unknown')

605     WRITE(6,*) '         '
        WRITE(6,'(A)') '$ despiking the individual mV : '
        WRITE(6,'(A)') '$            sigma cut at +/- : [3.0] '
	read(5,'(A)') sigchar

	if (sigchar(1:1) .eq. ' ') then
          sigma=3.0
	else
          read(sigchar,*) sigma
	end if

	write(50,*)' despiking data at +/-',sigma,'sigma level'	

        SIGMA = SIGMA * BBB     ! SETS UP CUT AT ? SIGMA FROM MEAN
        PLUS = AAA + SIGMA
        MINUS = AAA - SIGMA

        print*,' '
	print*,'clipping data outside ', plus,' to ',minus
	print*,'(shown as -...-...- line) '

        call pgslw(3)
	call pgsls(5)
	call pgmove(0.,plus)
	call pgdraw(real(oldnobs)+1.0,plus)
	call pgmove(0.,minus)
	call pgdraw(real(oldnobs)+1.0,minus)
	call pgslw(1)

	write(6,'(A)') '$ ... try another clipping level (Y/N) ? : [N] '
	read(5,'(A)') ans
	if (ans(1:1) .eq. 'Y' .or. ans(1:1) .eq. 'y')  goto 605
	call pgend

	DO I=1,EEE
             IF (CCC(I) .LT. PLUS .AND. CCC(I) .GT. MINUS) THEN
                 FFF = FFF + 1            ! COUNTS NO. OF POINTS IN NEW SAMPLE
                 WRITE(16,*) I, CCC(I), DDD(I)
             END IF
        END DO

           WRITE(6,*) '    '
           WRITE(6,*) EEE-FFF,' pairs removed from sample '
           WRITE(50,*) EEE-FFF,' pairs removed from sample '
           WRITE(6,*) '    '

	EEE=FFF 	      ! adjust no. of data to account for despiking

        CLOSE (UNIT=16)
        END


!-----------------------------------------------------------------------------!
!                          SUBROUTINE MANDESPIKE                              !
!-----------------------------------------------------------------------------!
 
!       DESPIKES ALL DATA SELECTED BY CURSOR 

       SUBROUTINE MANDESPIKE(CCC,DDD,EEE,DELPAIR,NODEL,SAMPLESIZE,option)

! CCC     REAL ARRAY(3000)    INDIVIDUAL DATA
! DDD     REAL ARRAY(3000)    INDIVIDUAL ERRORS
! EEE     INTEGER             NO OF DATA IN ORIGINAL SAMPLE
! DELPAIR INTEGER ARRAY (100) BAD DATA RETURNED FROM CURSOR
! NODEL   INTEGER             NO OF DELETED DATA
! SAMPLESIZE INTEGER          NO.OF DATA AFTER DESPIKING

	REAL CCC(3000), DDD(3000) 

        INTEGER I,EEE,NN,A,DELPAIR(300),NODEL,J,F,k, SAMPLESIZE, option

        NN = 1
        OPEN (UNIT=16,FILE='despiked.dat',STATUS='unknown')

c   sort deleted pair no.s into increasing order using straight 
c   insertion

        DO J=2,NODEL            !  picks out each element in turn
          A = DELPAIR(J)
          DO F = J-1,1,-1       !  look for the place to insert it 
            IF(DELPAIR(F).LE.A) GOTO 10
            DELPAIR(F+1)=DELPAIR(F)
          END DO
          F=0
10        DELPAIR(F+1)=A        !  insert it
        END DO
        
	DO I=1,EEE
             IF (I .NE. DELPAIR(NN)) THEN
                 WRITE(16,*) I, CCC(I), DDD(I)
             ELSE
                 NN = NN + 1             ! CYCLES THROUGH BAD DATA 
            ENDIF	  
        END DO

	EEE=EEE-NODEL

        CLOSE(UNIT=14)
        CLOSE(UNIT=16)

        END

!--------------------------------------------------------------------!
!                        SUBROUTINE KSTEST                           !
!--------------------------------------------------------------------!

        SUBROUTINE KSTEST(PAIR,MVOLTS_ZERO,ERR_ZERO,DATAFILE,NOBS,option)

c              KOLMOGOROV-SMIRNOV  NON-PARAMETRIC 2-SAMPLE TEST
c                      USING NAG ROUTINE  G08CDF

        implicit none

        INTEGER MAX_PAIRS               ! Size of array
        PARAMETER (MAX_PAIRS = 3000)

         real*4 pair(MAX_PAIRS),mvolts_zero(MAX_PAIRS),err_zero(MAX_PAIRS)
         real*4 subsample(500,500),suberr(500,500)
       
         REAL*4 SX(MAX_PAIRS),SY(MAX_PAIRS)
         REAL*4 X(MAX_PAIRS),Y(MAX_PAIRS), erry(MAX_PAIRS), D, P
         real*4 yinc(MAX_PAIRS),coadd(MAX_PAIRS),errcoadd(MAX_PAIRS),pset

         INTEGER IFAIL,N,M,NTYPE,size,nocoadd,dummy,total
         integer i,j,k,si,loopmin,loopmax,jsave,index,ii
         integer lastsize,datafile,saveindex,yy, nobs
         integer count, option

997    format(2x,a9,4x,a7,4x,a5,4x,a5,4x,a9,4x,a14)
998    FORMAT(3x,i4,9x,i4,5x,i4,7x,i2,6x,F8.4,6x,F8.4,6x,a6)

	if (datafile .eq. 14) then
           open(unit=14,file='mv_taucorrected.dat',status='old')        
	else
           open(unit=16,file='despiked.dat',status='old')        
	end if

	 write(6,*) ' '
	 write(6,*) ' ... running K-S test to check for statistical consistency'
	 write(6,*) '     in the  data. You are asked to choose the size of'
	 write(6,*) '     the subsamples (into which the data stream will be'
         write(6,*) '     cut) and the significance level at which pairs will'
	 write(6,*) '     be rejected and thrown out from the single scan or '
         write(6,*) '     concatonated data. '
993	 write(6,*) '  '
         write(6,'(A)') '$ subsample size  : '
         read(5,*) size 
	 if (size .ge. nobs) then
	   write(6,*) ' '
	   write(6,*) ' ... try again...'
	   write(6,*)' subsample size must be less than',nobs 
           goto 993
	 end if
         write(6,*) '  '
         write(6,'(A)') '$ significance level to reject data : '
         read(5,*) pset

	 write(50,*) ' KS-test :     subsample size =',size
	 write(50,*) '         : significance level =',pset

         jsave=0
         index=0
	          
         do j=1,100                   !  cycles through samples
             k=j-1
             loopmin=(size*k)+1
             loopmax=size*(k+1)
           do i=loopmin,loopmax       ! (1-30), (31-60), (61-90)
             read(datafile,*,end=600)pair(i),mvolts_zero(i),err_zero(i)
             index = index + 1        !  (1 - 30) local counter
	     saveindex = index
             subsample(j,index) = mvolts_zero(i)
             suberr(j,index) = err_zero(i)
          end do
             index=0
             jsave = jsave + 1         ! no. of (subsamples-1)

         end do

600         if (saveindex .le. size) then  !  watches for smaller  
            lastsize = saveindex           !  last subsample
            end if

            jsave = jsave + 1          ! no.of subsamples
          
!-----------------------------------------------------------------------

          write(6,*)' '
          write(6,997) 'subsample','coadded','total','ifail',
     :                'statistic','critical level'

          j=0
          nocoadd=0
          count=size		! counts total no. of data examined

1001       j=j+1                   ! follows the subsample (1,2,3,...)
            k=j-1
             loopmin=(size*k)+1
             loopmax=size*(k+1)
             index=0

          do dummy=loopmin,loopmax     ! (1-30), (31-60), (61-90)
             index = index + 1          !  (1 - 30)
             nocoadd = nocoadd + 1      ! total size of coadded array

             coadd(nocoadd) = subsample(j,index)  ! adds to data array
             errcoadd(nocoadd) = suberr(j,index)  ! adds to error array
             x(nocoadd)=coadd(nocoadd)       ! writes x array for NAG
           end do

	total=nocoadd

1002	   if (j+1 .eq. jsave) then
           size = lastsize
           end if

           do i=1,size
              y(i)=subsample(j+1,i)
	      erry(i)=suberr(j+1,i)
	      count=count+1		! counts total no. of data examined
           end do



!------------------------------------------------------------------
C                  PERFORM KS-TEST: Numerical recipes
!------------------------------------------------------------------

       IFAIL =  0     !

       n=nocoadd      ! no. in coadded sample
       m=size         ! no. in subsample

          CALL KSTWO(N, X, M, Y, D, P, SX, SY, IFAIL)

!------------------------------------------------------------------

c if p>0.2 then accept H(0), ie samples consistent, then concatonate

          if (p .gt. pset) then      ! passed  at percentage level p !
	     total=total+size     ! total size of coadded array          
             write(6,998) m,n,total,ifail,d,p,'passed'

	     do i=1,size
	     coadd(nocoadd+i)=y(i)   ! adds last subsample to coadded data
	     errcoadd(nocoadd+i)=erry(i)
	     end do

             if (count .eq. nobs) then
 		goto  1003
	     end if
             goto 1001
          else                                   ! failed !
  	     total=total     ! total size of coadded array
             write(6,998) m,n,total,ifail,d,p,'failed'
             j=j+1                  ! goto next subsample and re-try

             if (count .eq. nobs) then
 		goto  1003
	     end if
             goto 1002
          end if

!-----------------------------------------------------------------!
C   CLEAR UP and OUTPUT DATA BACK TO SUBROUTINE 'MEAN', 'STERR'   !
!-----------------------------------------------------------------!
 
1003   open (unit=20,file='coadded.dat',status='unknown')

       do i=1,total
         write(20,*) i, coadd(i), errcoadd(i)
       end do

       nobs=total
       datafile=20

       close (unit=14)
       close (unit=16)
       close (unit=20)
   
       END

*************************************************************
***   Global Section Definition File (GSDF) Subroutines   ***
*************************************************************

      SUBROUTINE GET_C13DAT (FILENAME, QARRAY, RARRAY, DCOUNT, STATUS)
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*  Import:
      CHARACTER*64 FILENAME
      INTEGER KROWCOL
*  Import/Export:
      INTEGER DCOUNT
*  Export:
      real qarray(*)
      real rarray(*)
*  status:
      integer status
*  Local variables:
      character*16 dname
      real*8 d1,d2,d3
      real*4 r1,r2
      real cellx,celly
      integer starts(gsd__mxdim)
      integer ends(gsd__mxdim)
      integer udimvals(gsd__mxdim)
      integer uactdims
      integer i,k,kk,kcount
* 
      if (status.ne.adam__ok) return
      if (dcount.le.0) return
*
*  Get cell X size
*
      dname='c6dx'
      status=adam__ok
      call get_double (filename, dname, cellx, status)
*
*  Read Pointing History array to determine whether X fit or Y fit.
*
*  First, inquire file for dimension values (UDIMVALS) of c14phist array.
*  This is indicated by setting starts and ends to 0.

      dname='c14phist'
      kcount = 0
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, qarray, kcount, status)
      starts(1) = 1
      starts(2) = 1
      ends(1)   = udimvals(1)
      ends(2)   = udimvals(2)
*
* Now, read an entire row or column of MAP data
*
      kcount = 2*dcount
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, qarray, kcount, status)
*
* Now, read raster or grid data
* First, get dimension values (UDIMVALS) of the C13DAT array
*
      dname='c13dat'
      kcount = dcount
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, rarray, kcount, status)

      starts(1) = 1
      starts(2) = 1
      ends(1)   = udimvals(1)
      ends(2)   = udimvals(2)
      kcount    = dcount
*
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, rarray, kcount, status)
*
*  Now, reconstruct X and Y data.
*  KCOUNT is actual amount of data returned in RARRAY
*
      k = 0
      do kk = 1, 2*kcount-1, 2
         k = k+1
         qarray(k) = qarray(kk) * cellx
      enddo
*
      dcount = kcount
*
999   continue
*
      end
* 

      SUBROUTINE GET_C13SPV (FILENAME, QARRAY, DARRAY, DCOUNT, STATUS)
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*  Import:
      character*64 filename
      integer krowcol
*  Import/Export:
      integer dcount
*  Export:
      real qarray(*)
      real*8 darray(*)
*  Status:
      integer status
*  Local variables:
      character*16 dname
      real*8 d1,d2,d3
      real*4 r1,r2
      real cellx,celly
      integer starts(gsd__mxdim)
      integer ends(gsd__mxdim)
      integer udimvals(gsd__mxdim)
      integer uactdims
      integer i,k,kk,kcount
* 
      if (status.ne.adam__ok) return
      if (dcount.le.0) return

*
*  Get cell X size
*
      dname='c6dx'
      status=adam__ok
      call get_double (filename, dname, cellx, status)
*
*  Read Pointing History array to determine whether X fit or Y fit.
*
*  First, inquire file for dimension values (UDIMVALS) of c14phist array.
*  This is indicated by setting starts and ends to 0.

      dname='c14phist'
      kcount = 0
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, qarray, kcount, status)
      starts(1) = 1
      starts(2) = 1
      ends(1)   = udimvals(1)
      ends(2)   = udimvals(2)
*
* Now, read an entire row or column of MAP data
*
      kcount = 2*dcount
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, qarray, kcount, status)
*
* Now, read raster or grid data
* First, get dimension values (UDIMVALS) of the C13DAT array
*
      dname='c13spv'
      kcount = dcount
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, darray, kcount, status)

      do i=1,uactdims
         starts(i)=1
         ends(i)= udimvals(i)
      enddo

      kcount    = dcount
*
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, darray, kcount, status)
*
*  Now, reconstruct X and Y data.
*  KCOUNT is actual amount of data returned in DARRAY
*
      k = 0
      do kk = 1, 2*kcount-1, 2
         k = k+1
         qarray(k) = qarray(kk) * cellx
      enddo
*
      dcount = kcount
*
999   CONTINUE
*
      END
* 
*

      SUBROUTINE GET_C13RAW_ERROR (FILENAME, QARRAY, DARRAY, DCOUNT, STATUS)
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*  Import:
      character*64 filename
      integer krowcol
*  Import/Export:
      integer dcount
*  Export:
      real qarray(*)
      real*8 darray(*)
*  Status:
      integer status
*  Local variables:
      character*16 dname
      real*8 d1,d2,d3
      real*4 r1,r2
      real cellx,celly
      integer starts(gsd__mxdim)
      integer ends(gsd__mxdim)
      integer udimvals(gsd__mxdim)
      integer uactdims
      integer i,k,kk,kcount
* 
      if (status.ne.adam__ok) return
      if (dcount.le.0) return
*
*  Get cell X size
*
      dname='c6dx'
      status=adam__ok
      call get_double (filename, dname, cellx, status)
*
*  Read Pointing History array to determine whether X fit or Y fit.
*
*  First, inquire file for dimension values (UDIMVALS) of c14phist array.
*  This is indicated by setting starts and ends to 0.

      dname='c14phist'
      kcount = 0
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, qarray, kcount, status)
      starts(1) = 1
      starts(2) = 1
      ends(1)   = udimvals(1)
      ends(2)   = udimvals(2)
*
* Now, read an entire row or column of MAP data
*
      kcount = 2*dcount
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, qarray, kcount, status)
*
* Now, read raster or grid data
* First, get dimension values (UDIMVALS) of the C13DAT array
*
      dname='c13raw_error'
      kcount = dcount
      status=adam__ok
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, 0, 0, darray, kcount, status)

      do i=1,uactdims
         starts(i)=1
         ends(i)= udimvals(i)
      enddo

      kcount    = dcount
*
      call get_gsd_data (filename, dname,
     :     udimvals, uactdims, starts, ends, darray, kcount, status)
*
*  Now, reconstruct X and Y data.
*  KCOUNT is actual amount of data returned in DARRAY
*
      k = 0
      do kk = 1, 2*kcount-1, 2
         k = k+1
         qarray(k) = qarray(kk) * cellx
      enddo
*
      dcount = kcount
*
999   CONTINUE
*
      END
* 
*


*
      SUBROUTINE GET_CHARACTER (FILENAME, CNAME, CVALUE, STATUS)
*  Get a character item from GSD file
*  The CVALUE is a Character*16 data item.
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*  Define GET_GSD_DATA arguments
*  Imports:
      CHARACTER*64 FILENAME
      CHARACTER*16 CNAME
*  Export:
*      CHARACTER*16 CVALUE
      CHARACTER*16 CVALUE
*  Status:
      INTEGER STATUS
*  Local:
      INTEGER DUMMY
      INTEGER DDD
      INTEGER DD
      INTEGER DDDB
      INTEGER DDDD
* 
      INTEGER NUMBER
      INTEGER FD
      REAL VERSION
      CHARACTER*1 TYPE
      character*(gsd__szunit) unit
      character*20 label        !gsd file label
      integer nitem	
      logical array
      integer index (gsd__szindex)	
      IF (STATUS.NE.ADAM__OK) RETURN

C      CALL GET_GSD_DATA (FILENAME, CNAME,
C     :     DUMMY, DDDB, DDDD, DD, CVALUE, DDD, STATUS)
      call gsd_open_read (filename, fd, version, label, nitem, status)
      CALL GSD_FIND (fd, cname, number, unit, type,
     :     array, index, status)
      IF (status .ne. adam__ok) return
      call gsd_get0c (index,cvalue, status)

      CALL GSD_CLOSE(FD)

      END
* 
      SUBROUTINE GET_REAL (FILENAME, RNAME, RVALUE, STATUS)
*  Get Real valued scalar from GSD file
*  The RVALUE is Real*4 data item.
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
*  Define GET_GSD_DATA arguments
*  Imports:
      CHARACTER*64 FILENAME
      CHARACTER*16 RNAME
*  Export:
      REAL*4 RVALUE
*  Status:
      INTEGER STATUS
*  Local:
      INTEGER DUMMY
* 
      IF (STATUS.NE.ADAM__OK) RETURN
*
      CALL GET_GSD_DATA (FILENAME, RNAME,
     :     DUMMY, DUMMY, DUMMY, DUMMY, RVALUE, DUMMY, STATUS)
      END
*
      SUBROUTINE GET_INTEGER (FILENAME, INAME, IVALUE, STATUS)
*  Get Integer*4 valued scalar from GSD file
*  The IVALUE is Integer*4 data item.
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
*  Define GET_GSD_DATA arguments
*  Imports:
      CHARACTER*64 FILENAME
      CHARACTER*16 INAME
*  Export:
      INTEGER IVALUE 
*  Status:
      INTEGER STATUS
*  Local:
      INTEGER DUMMY
* 
      IF (STATUS.NE.ADAM__OK) RETURN
*
      CALL GET_GSD_DATA (FILENAME, INAME,
     :     DUMMY, DUMMY, DUMMY, DUMMY, IVALUE, DUMMY, STATUS)
      END
*
      SUBROUTINE GET_LOGICAL (FILENAME, LNAME, LVALUE, STATUS)
*  Get Integer*4 valued scalar from GSD file
*  The IVALUE is Integer*4 data item.
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
*  Define GET_GSD_DATA arguments
*  Imports:
      CHARACTER*64 FILENAME
      CHARACTER*16 LNAME
*  Export:
      LOGICAL LVALUE 
*  Status:
      INTEGER STATUS
*  Local:
      INTEGER DUMMY
* 
      IF (STATUS.NE.ADAM__OK) RETURN
*
      CALL GET_GSD_DATA (FILENAME, LNAME,
     :     DUMMY, DUMMY, DUMMY, DUMMY, LVALUE, DUMMY, STATUS)
      END
*
*
      SUBROUTINE GET_DOUBLE (FILENAME, DNAME, DVALUE, STATUS)
*  Get Real*8 valued scalar from GSD file
*  The DVALUE is Real*8 data item.
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
*  Define GET_GSD_DATA arguments
*  Imports:
      CHARACTER*64 FILENAME
      CHARACTER*16 DNAME
*  Export:
      REAL*8 DVALUE
*  Status:
      INTEGER STATUS
*  Local:
      INTEGER DUMMY
* 
      IF (STATUS.NE.ADAM__OK) RETURN
*
      CALL GET_GSD_DATA (FILENAME, DNAME,
     :     DUMMY, DUMMY, DUMMY, DUMMY, DVALUE, DUMMY, STATUS)
      END
*
*
      SUBROUTINE GET_AIRMASS (FILENAME, AIRMASS, STATUS)
*+ Get Airmass from GSD file
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*  Import:
      CHARACTER*64 FILENAME
*  Export:
      REAL*4 AIRMASS
*  Status:
      INTEGER STATUS
*  Local variables:
      REAL*4 C12SCAN_TABLE_1(6)
      CHARACTER*16 DNAME
      INTEGER STARTS(GSD__MXDIM)
      INTEGER ENDS(GSD__MXDIM)
      INTEGER UDIMVALS(GSD__MXDIM)
      INTEGER UACTDIMS
      INTEGER KCOUNT
* 

      IF (STATUS.NE.ADAM__OK) RETURN

      DNAME='C12SCAN_TABLE_1'
      STARTS(1)=1
      STARTS(2)=1
      ENDS(1)=1
      ENDS(2)=6
      KCOUNT=6
      STATUS=ADAM__OK
*  Note: D1 is Real*8
      CALL GET_GSD_DATA (FILENAME, DNAME,
     :     UDIMVALS, UACTDIMS, STARTS, ENDS, C12SCAN_TABLE_1, KCOUNT, STATUS)
*
      AIRMASS = C12SCAN_TABLE_1(2)
*
      END
*
*+  GET_GSD_DATA
      SUBROUTINE GET_GSD_DATA (FILENAME, DNAME,
     :        UDIMVALS, UACTDIMS, STARTS, ENDS,
     :        RVALUES, DCOUNT, STATUS)
*    Description :
*        Read GSD file <FILENAME>, search for item <DNAME>, and tranfer
*      <DCOUNT> amount of data into RVALUES.
*      Based on routines developed by Jon Fairclough.
*
*        If <DNAME> is an array, its dimensional range is returned
*      in <UDIMVALS> and the actual number of dimensions is returned
*      in<UACTDIMS>. <STARTS> and <ENDS> are dimensioned variables which
*      indicate the starting and ending element numbers to transfer into
*      array <RVALUES>.
*
*      Note: RVALUES is really a pointer to storage for any item requested.
*            It could point to items of type real, double precision,
*            or character, as well as arrays of these types.
*
*        If <DNAME> is a scalar, acceptable types are real*4 (R), real*8 (D),
*      integer*4 (I) and character*(*) (C).
*        
*    Authors :
*      Firmin J. Oliveira
*    History :
*      12-Mar-1991 : Original
*    Type Definitions :
      IMPLICIT NONE
*    Global constants :
       include '/star/include/adam_err'
       include '/star/include/prm_par'
       include 'gsd_pars.inc'
*    Import :
      character*(*) filename
      character*(*) dname
      integer starts(*)
      integer ends(*)
*    Export :
      integer udimvals(*)
      integer uactdims
      real rvalues (*)
      integer dcount
*    Status :
      integer status 
*    Local variables :
      integer i			! the ususal running variable
* GSD_OPEN_READ
      integer fd	        !number of gsd file
      real version	        !gsd file version number
      character*20 label        !gsd file label
      integer nitem	        ! number of gsd items in file
* GSD_FIND
      character*(gsd__szname) name	!name of gsd item sought
      integer number			!number of gsd item
      character*(gsd__szunit) unit 	!unit of gsd item
      character*1 type			!type of gsd item
      logical array			!true if item array
      integer index (gsd__szindex)	!index for accessing gsd item
* GSD_INQ_SIZE (Replaces GSD_INQ_ARRAY)
      integer maxdims                   ! maximum number of dimensions
      integer dimvals (gsd__mxdim) 	! values of dim. variables
c     INTEGER DIMNUMBERS(GSD__MXDIM)     !The ordinal numbers of dim. scalars
      character*20 dimnames (gsd__mxdim) ! names of dimensioning variables
      character*20 dimunits (gsd__mxdim) ! units of dim. variables
c      INTEGER DIMINDEX(GSD__SZINDEX, GSD__MXDIM) !indices of dim. quantities
      integer actvals                   ! actual number of values
      integer actdims			! actual number of dimensions
      integer size			! total size in cells
* GSD_GET1R
      integer actvals             !count of actual data elems. transfered
* Special handling
      character*64 old_filename
      save old_filename
      data old_filename/' '/
*
*-
*     
      IF (STATUS .NE. ADAM__OK) RETURN

* open GSD file (if it is not open):

      IF (filename.ne.old_filename) then
         call gsd_open_read (filename, fd, version, label, nitem,status)
         if (status .ne. adam__ok) return
         old_filename=filename
      ENDIF

* find the relevant items:

      CALL GSD_FIND (fd, dname, number, unit, type,
     :             array, index, status)
      IF (status .ne. adam__ok) return

      IF (array) then
*  Get dimensional info on data array
         call gsd_inq_size(fd,number,gsd__mxdim,dimnames,dimunits,
     :    dimvals,actdims,size,status)

c         CALL GSD_INQUIRE_ARRAY(INDEX, MAXDIMS,
c     :             DIMNUMBERS, DIMNAMES, DIMUNITS,
c     :             DIMVALS, DIMINDEX, ACTDIMS, SIZE, STATUS)

*  If ENDS(1)=0 then transfer DIMVALS and ACTDIMS only
         IF (ends(1).le.0) then
            DO i=1,actdims
               udimvals(i) = dimvals(i)
            ENDDO
            uactdims = actdims
            dcount=0
            return            
         ENDIF

         IF (type.eq.'R') then
            call gsd_get1r (index, actdims, dimvals, starts, ends, 
     :             rvalues, actvals, status)
         ELSEIF (type.eq.'D') then
            call gsd_get1d (index, actdims, dimvals, starts, ends,
     :             rvalues, actvals, status)

         ELSE
* Currently, unsupported type
            status = -1
            actvals = 0
         ENDIF

*  Return total count of actual data elements transfered
         dcount = actvals

      ELSE

* GSD returns type as an upper case character...

*  Get a scalar item
         IF (type.eq.'R') then
            call gsd_get0r (index, rvalues(1), status)
            dcount = 1
         ELSEIF (type.eq.'D') then
            call gsd_get0d (index, rvalues(1), status)
            dcount = 1
         ELSEIF (type.eq.'C') then     
            call gsd_get0c (index, rvalues(1), status)
            dcount = 1
         ELSEIF (type.eq.'I') then
            call gsd_get0i (index, rvalues(1), status)
            dcount = 1
         ELSEIF (type.eq.'L') then
            call gsd_get0l (index, rvalues(1), status)
            dcount = 1
         ELSE
            status = -1
         ENDIF
      ENDIF

      END      

*-----------------------------------------------------------------------

      subroutine uucase (string)

*  Subroutine to convert entire string to upper case using available
*  system services.

      implicit    none

*     Formal parameters:

      character*(*)  string

*     Local variables:

      integer     i
      integer     ich

*  Ok, go...

c     type *, ' --- uucase ---'
c     type *, '     input string  = ', string

      do i = 1, len (string)
        ich = ichar (string(i:i))
        if (ich.ge.097 .and. ich.le.122) string(i:i) = char (ich-32)
      end do
 
C     TYPE *, '     output string = ', string


      return
      end

*-----------------------------------------------------------------------

      subroutine ulcase (string)

*  Subroutine to convert entire string to lower case using available
*  system services.

      implicit    none

*     Formal parameters:

      character*(*)  string

*     Local variables:

      integer     i
      integer     ich

*  Ok, go...

c     type *, ' --- ulcase ---'
c     type *, '     input string  = ', string

      do i = 1, len (string)
        ich = ichar (string(i:i))
        if (ich.ge.065 .and. ich.le.90) string(i:i) = char (ich+32)
      end do
 
c     type *, '     output string = ', string


      return
      end


******** Now include the NR KS test from the SCUBA software

********* These are the numerical recipes routines ******************
 
* This is the two sided KS statistic
 
      SUBROUTINE KSTWO(N1, DATA1, N2, DATA2, D, PROB, SORT1, SORT2, 
     :     STATUS)
*+
*  Name:
*     KSTWO
 
*  Purpose:
*     Returns the probability that two data sets are from the same sample
 
*  Language:
*     Starlink Fortran 77
 
*  Invocation:
*     CALL KSTWO( N1, DATA1, N2, DATA2, D, PROB, SORT1, SORT2, STATUS )
 
*  Description:
*     This routine computes the two sided KS statistic of two arrays
*     and returns the probability that the two arrays are drawn from the
*     same sample.
 
*  Arguments:
*     N1 = INTEGER (Given)
*        Number of elements in first array
*     DATA1 ( N1 ) = REAL (Given & Returned)
*        First input array
*     N2 = INTEGER (Given)
*        Number of elements in second array
*     DATA2 ( N2 ) = REAL (Given & Returned)
*        Second input array
*     D  = REAL (Returned)
*        Maximum distance between cumulative distribution functions
*     PROB = REAL (Returned)
*        Probability two arrays are from the same sample
*     SORT1 = REAL (Returned)
*        Sorted form of array DATA1
*     SORT2 = REAL (Returned)
*        Sorted form of array DATA2
*     STATUS = INTEGER (Given & Returned)
*        Global status value
 
*  Notes:
 
*  Algorithm:
*     This routine first sorts both arrays and then compares cumulative
*     distribution functions. The largest separation between these functions
*     is then used to calculate the Kolmogorov-Smirnov statistic.
 
*  References:
*     - Press et al, 1992, "Numerical Recipes in FORTRAN", 2nd edition (CUP)
 
*  Implementation status:
*     - Bad pixels must be removed before passing to this routine
 
*  Authors:
*     TIMJ: Tim Jenness (JACH)
*     {enter_new_authors_here}
 
*  History:
*     1996 October 22 (TIMJ):
*       Original Starlink version
*     {enter_changes_here}
 
*  Bugs:
*     {note_any_bugs_here}
 
*-
 
*  Type Definitions:
      IMPLICIT  NONE
 
*  Global Constants:
      INCLUDE '/star/include/sae_par'               ! SSE global definitions
 
*  Arguments Given:
      INTEGER N1
      INTEGER N2
 
*  Arguments Given and Returned:
      REAL    DATA1( N1 )
      REAL    DATA2( N2 )
 
*  Arguments Returned:
      REAL    D
      REAL    PROB
      REAL    SORT1( N1 )
      REAL    SORT2( N2 )
 
*  Status:
      INTEGER STATUS             ! Global status
 
*  External References:
      EXTERNAL PROBKS
      REAL     PROBKS                 ! KS statistic
 
*  Local Variables:
      REAL    D1                      ! Value from first set
      REAL    D2                      ! Value from second set
      REAL    DTEMP                   ! Difference between FN2 and FN1
      REAL    EN                      ! Weighted number of points
      REAL    FN1                     ! Fraction through data set 1
      REAL    FN2                     ! Fraction through data set 2
      INTEGER I                       ! Loop counter
      INTEGER J1                      ! Counter
      INTEGER J2                      ! Counter
      REAL    R1                      ! REAL(N1)
      REAL    R2                      ! REAL(N2)
*.
 
*     Check the global inherited status.
      IF (STATUS .NE. SAI__OK) RETURN 
 
*     Copy data to the scratch arrays
      DO I = 1, N1
         SORT1(I) = DATA1(I)
      END DO
      DO I = 1, N2
         SORT2(I) = DATA2(I)
      END DO
 
 
*     Sort the data before further processing
 
      CALL KPG1_QSRTR(N1, 1, N1, SORT1, STATUS)
      CALL KPG1_QSRTR(N2, 1, N2, SORT2, STATUS)
 
      IF (STATUS .NE. SAI__OK) THEN 
         PROB = -1.0
         D    = 0.0
         RETURN 
      END IF
 
*     Initialise variables
      R1 = REAL(N1)
      R2 = REAL(N2)
 
      J1=1
      J2=1
      FN1=0.
      FN2=0.
      D=0.
 
*     Loop through data
      DO WHILE (J1.LE.N1 .AND. J2.LE.N2)
 
         D1 = SORT1(J1)
         D2 = SORT2(J2)
 
         IF (D1 .LE. D2) THEN
            FN1 = J1 / R1
            J1 = J1+1
         END IF
 
         IF (D2 .LE. D1) THEN
            FN2 = J2 / R2
            J2 = J2 + 1
         END IF
 
*        Find distance between the points
         DTEMP = ABS(FN2 - FN1)
         IF (DTEMP.GT.D) D = DTEMP
      END DO
 
*     Find KS statistic associated with maximum separation
      EN = SQRT( R1 * R2 / (R1 + R2) )
      PROB = PROBKS(( EN + 0.12 + 0.11 / EN ) * D)
 
      END
 
 
* This is the Kolmogorov-Smirnov probability distribution
 
      REAL FUNCTION PROBKS( ALAM )
*+
*  Name:
*     PROBKS
 
*  Purpose:
*     Returns the KS statistic
 
*  Language:
*     Starlink Fortran 77
 
*  Invocation:
*     PROB = PROBKS( ALAM )
 
*  Arguments:
*     PROBKS = REAL (Returned via FUNCTION)
*       The probability the samples were identical
*     ALAM   = REAL (Given)
*       Measure of separation of two cumulative distributions
 
*  References:
*     - Press et al, 1992, "Numerical Recipes in FORTRAN", 2nd edition (CUP)
 
 
*  Authors:
*     TIMJ: Tim Jenness (JACH)
*     {enter_new_authors_here}
 
*  History:
*     1996 October 22 (TIMJ):
*       Original Starlink version
*     {enter_changes_here}
 
*  Bugs:
*     {note_any_bugs_here}
 
*-
 
*  Type Definitions:
      IMPLICIT NONE
 
*  Arguments given:
      REAL    ALAM                  ! Dimensionless distance
 
*  Local constants:
      REAL    EPS1                  ! Required precision
      PARAMETER (EPS1 = 0.001)
      REAL    EPS2                  !
      PARAMETER (EPS2 = 1.E-8)
 
*  Local variables:
      REAL    A2                    !
      REAL    FAC                   !
      INTEGER J                     ! Loop counter
      REAL    TERM                  !
      REAL    TERMBF                !
*.
 
*  Initialise loop variables
      A2 = -2. * DBLE(ALAM)**2
      FAC = 2.
      PROBKS = 0.0
      TERMBF = 0.0
 
*  Begin loop
      DO J = 1, 100
 
         TERM = FAC * EXP( A2 * J**2 )
         PROBKS = PROBKS + TERM
 
*     Return if term is smaller than precision
         IF (ABS(TERM) .LE. EPS1 * TERMBF .OR. 
     :        ABS(TERM).LE. EPS2 * PROBKS) THEN
            RETURN
         END IF
 
        FAC = -FAC
        TERMBF = ABS(TERM)
 
      END DO
 
      PROBKS=1.
      RETURN
      END


****** Sort routine from KAPPA

      SUBROUTINE KPG1_QSRTR( EL, LOW, HIGH, ARRAY, STATUS )
*+
*  Name:
*     KPG1_QSRTX
 
*  Purpose:
*     Sorts a vector via the Quicksort algorithm.
 
*  Language:
*     Starlink Fortran 77
 
*  Invocation:
*     CALL KPG1_QSRTx( EL, LOW, HIGH, ARRAY, STATUS )
 
*  Description:
*     This routine sorts a vector in situ between an upper and lower
*     bounds using the Quicksort algorithm.
 
*  Arguments:
*     EL = INTEGER (Given)
*        The number of elements in the array that is to be sorted.
*     LOW = INTEGER (Given)
*        The lower bound within the array, below which the array
*        elements will not be sorted.  It should be less than the
*        upper bound and must be within the array.  In the latter case
*        an error will result and the routine will not sort the array.
*     HIGH = INTEGER (Given)
*        The upper bound within the array, above which the array
*        elements will not be sorted.  It should be greater than the
*        lower bound and must be within the array.  In the latter case
*        an error will result and the routine will not sort the array.
*     ARRAY( EL ) = ? (Given and Returned)
*        The array to be sorted.
*     STATUS = INTEGER (Given and Returned)
*        The global status.
 
*  Algorithm:
*     Quicksort works by picking a random "pivot" element in the array
*     then moving every element that is bigger to one side of the
*     pivot, and every element that is smaller to the other side.  The
*     procedure is repeated with the two subdivisions created by the
*     pivot.  When the number of elements in a subdivision reaches two,
*     the array is sorted.
*
*     Since recursion is not possible in Fortran, pushdown stacks are
*     used to mimic the recursive operation of the Quicksort algorithm.
*     These are also more efficient than recursion because they avoid
*     the expensive subroutine calls, especially for the many small
*     subdivisions.  The stacks contains the subdivisions to be sorted.
*     The stack is popped to obtain a subfile to sort.  The partitioning
*     pushs the larger subdivisions on to the stack, and the smaller
*     subdivision is processed immediately.  Hence the stack size
*     is only lg(EL).
 
*  Implementation Status:
*     -  There is a routine for each of the data types integer, real,
*     double precision, and character: replace "x" in the routine nam
*     by I, R, D, or C respectively as appropriate.
*     -  If the maximum bound is less than the minimum, the bounds are
*     swapped.
 
*  References:
*     -  Sedgwick, R., 1988, "Algorithms" (Addison-Wesley).
 
*  Timing:
*     For N elements to be sorted the timing goes as NlnN.
 
*  Authors:
*     MJC: Malcolm J. Currie (STARLINK)
*     {enter_new_authors_here}
 
*  History:
*     1991 January 11 (MJC):
*        Original version.
*     {enter_changes_here}
 
*  Bugs:
*     {note_any_bugs_here}
 
*-
 
*  Type Definitions:
      IMPLICIT NONE              ! No implicit typing
 
*  Global Constants:
      INCLUDE '/star/include/sae_par'          ! Standard SAE constants
 
*  Arguments Given:
      INTEGER
     :  EL,
     :  LOW,
     :  HIGH
 
*  Arguments Given and Returned:
      REAL
     :  ARRAY( EL )
 
*  Status:
      INTEGER STATUS             ! Global status
 
*  Local Constants:
      INTEGER MXSTAK             ! The size of the stack which is
                                 ! log base 2 of the maximum number of
                                 ! elements in the array
      PARAMETER ( MXSTAK = 32 )
 
*  Local Variables:
      INTEGER
     :  I,                       ! Ascending pointer to an array element
     :  J,                       ! Descending pointer to an array
                                 ! element
     :  LOWER( MXSTAK ),         ! Stack for the elements of the array
                                 ! below the pivot element
     :  LBND,                    ! Polarity-checked version of the lower
                                 ! bound
     :  PSTACK,                  ! Pointer to the stack
     :  UBND,                    ! Polarity-checked version of the upper
                                 ! bound
     :  UPPER( MXSTAK )          ! Stack for the elements of the array
                                 ! above the pivot element
 
      REAL
     :  PIVOT,                   ! Pivot element
     :  TEMP                     ! Used for swapping array elements
 
*.
 
*    Check inherited global status.
 
      IF ( STATUS .NE. SAI__OK ) RETURN
 
*    Validate the array limits.
*    ==========================
 
*    Check that they lie within the array's bounds.  If not report an
*    the error and exit.
 
      IF ( LOW .LT. 1 .OR. LOW .GT. EL .OR. HIGH .LT. 1 .OR.
     :     HIGH .GT. EL ) THEN
         STATUS = SAI__ERROR
         print *, 'Sorting limits are outside bounds of array'
         GOTO 999
      END IF
 
*    Check the polarity.
 
      IF ( LOW .GT. HIGH ) THEN
 
*       A swap is necessary.
 
         LBND = LOW
         UBND = HIGH
      ELSE
 
*       Just copy the input values.
 
         LBND = LOW
         UBND = HIGH
      END IF
 
*    ^^^^^^^^^^^^^^^^^^^^^^^^^^
 
*    Intialise the stacks.
 
      LOWER( 1 ) = LBND
      UPPER( 1 ) = UBND
      PSTACK = 1
 
*    Loop until the stack is empty.
 
      DO WHILE ( PSTACK .GT. 0 )
 
*       Pop the stack.
 
         IF ( LOWER( PSTACK ) .GE. UPPER( PSTACK ) ) THEN
            PSTACK = PSTACK - 1
         ELSE
 
*          Partition the array.
*          ====================
 
            I = LOWER( PSTACK )
            J = UPPER( PSTACK )
            PIVOT = ARRAY( J )
 
*          Move in from both sides towards the pivot element.
 
            DO WHILE ( I .LT. J )
               DO WHILE ( ( I .LT. J )  .AND. ARRAY( I ) .LE. PIVOT )
                  I = I + 1
               END DO
 
               DO WHILE ( ( J .GT. I )  .AND. ARRAY( J ) .GE. PIVOT )
                  J = J - 1
               END DO
 
*             If the pivot element is not yet reached, it means that two
*             elements on either side are out of order, so swap them.
 
               IF ( I .LT. J ) THEN
                  TEMP = ARRAY( I )
                  ARRAY( I ) = ARRAY( J )
                  ARRAY( J ) = TEMP
               END IF
            END DO
 
*          Move the pivot element back to its proper place in the array.
 
            J = UPPER( PSTACK )
            TEMP = ARRAY( I )
            ARRAY( I ) = ARRAY( J )
            ARRAY( J ) = TEMP
 
*          Push values on to the stacks to further subdivide the
*          array.
 
            IF ( ( I - LOWER( PSTACK ) ) .LT.
     :           ( UPPER( PSTACK ) - I ) ) THEN
               LOWER( PSTACK + 1 ) = LOWER( PSTACK )
               UPPER( PSTACK + 1 ) = I - 1
               LOWER( PSTACK )     = I + 1
            ELSE
               LOWER( PSTACK + 1 ) = I + 1
               UPPER( PSTACK + 1 ) = UPPER( PSTACK )
               UPPER( PSTACK )     = I - 1
            END IF
 
*          Increment the stack counter.
 
            PSTACK = PSTACK + 1
         END IF
      END DO
 
  999 CONTINUE
 
      END
