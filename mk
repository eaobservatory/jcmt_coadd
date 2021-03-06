
#  N.B. the previous line should be blank.
#+
#  Name:
#     mk

#  Purpose:
#     Invoke make to build and install an application package.

#  Type of Module:
#     Shell script.

#  Description:
#     This script should normally be used to invoke the make utility to
#     build and install an application package and to perform other
#     housekeeping tasks. It invokes the make utility after first
#     defining appropriate environment variables and macros for the
#     computer system in use.

#  Invocation:
#     The user of this script should normally first define the SYSTEM
#     environment variable to identify the host computer system (see
#     the "Supported Systems" section). This script should then be used
#     in the same way as the make utility would be used. For instance,
#     to build, install and test the package, you might use the following
#     commands:
#
#        mk build
#	 mk install
#	 mk test
#	 mk clean
#
#     See the makefile prologue for further details.

#  Supported Systems:
#     The following UNIX systems are currently supported and may be
#     identified by defining the SYSTEM environment variable
#     appropriately before invoking this script:
#        alpha_OSF1
#           DEC Alpha machines running OSF1
#        mips
#           DECstations running Ultrix
#        sun4
#           SUN Sparcstations running SunOS (4.x)
#        sun4_Solaris
#           SUN Sparcstations running SunOS 5.x (Solaris)
#
#     This script will exit without action if the SYSTEM environment
#     variable is not defined. A warning will be issued if it is
#     defined but is not set to one of the values above. In this case,
#     no additional environment variables will be defined by this
#     script (any that are pre-defined will be passed on to make
#     unaltered).

#  Notes on Porting:
#     If your machine or system setup does not appear in this script,
#     then it should be possible to build and install the package by adding a
#     new case to this script with appropriate definitions (probably
#     based on one of the existing implementations).

#  make Macros:
#     The following "global" make macros are used in the associated
#     makefile and may be changed by means of appropriate environment
#     variable definitions (in each case the default is shown in
#     parentheses). Note that these macros are provided to allow
#     external control over the directories in which software is
#     installed, etc., so they should not normally be re-defined within
#     this script.
#     
#        STARLINK (/star)
#	    Pathname of the root directory beneath which Starlink
#	    software is currently installed. This indicates to the build
#           process where to find other Starlink software (include files,
#	    libraries, etc.) which it uses.
#        INSTALL ($HOME)
#	    Pathname of the root directory beneath which the package will be
#	    installed for use. Your home directory will be used by
#	    default. This macro is provided to allow the package to be
#	    installed locally for personal use (e.g. during development
#	    or testing). It should be set to the $STARLINK directory if
#	    you want to add the package into an already installed set of
#	    Starlink software. You should ensure that the appropriate
#	    sub-directories appear on any relevant search paths which
#	    your system uses for locating software (e.g. binaries and
#	    libraries).
#        EXPORT (.)
#	    Pathname of the directory into which compressed tar files
#	    will be written if the "export" or "export_source" make
#	    targets are used to produce an exportable copy of the package or
#	    its source files. The current working directory (i.e. the
#	    package's source directory) will be used by default.
#
#     The following "local" make macros are used in the associated
#     makefile and should normally be over-ridden (when appropriate) by
#     environment variable definitions within this script. In each case
#     the default is shown in parentheses.
#
#        ALINK (alink)
#           The command used to link an atask. Usually this is the 
#           alink command.
#           absent).
#        AR_IN (ar r)
#	    The command to use to insert an object (.o) file into an
#	    archive (.a) library. On some systems the variation 'ar -r'
#	    may be required instead.
#        CC (cc)
#	    The C compiler command to use. CCDPACK does currently
#	    require a C compiler.
#        CFLAGS (-O)
#	    The C compiler options to use.
#        FC (f77)
#	    The Fortran compiler command to use. CCDPACK requires a Fortran
#	    77 compiler that supports the common "permitted" Starlink
#	    extensions, as documented in Starlink General Paper SGP/16.
#	    (These include only the most common extensions, such as
#	    long names, end of line comments, include files, etc. together
#           with the %VAL mechanism)
#        FFLAGS (-O)
#	    The Fortan compiler options to be used.
#        LINK (ln)
#	    The command required to establish a link to a file. The
#	    default assumes POSIX.2 behavior, which only provides a
#	    "hard" link operating within a single file system. If the
#	    host operating system allows "symbolic" links, then this
#	    macro might be re-defined as 'ln -s'. Alternatively, if the
#	    use of multiple file systems is essential but not supported
#	    by any form of link, then a copy command could be
#	    substituted (e.g. 'cp -p'), at some cost in file space.
#        RANLIB (echo >/dev/null)
#	    The command required to "randomise" an object library. By
#	    default, this operation is not performed (the default acts
#	    as a null command). On systems which require it, this
#	    should typically be set to 'ranlib'.
#	 TAR_IN (pax -w -v -x ustar -f)
#	    Command to use to insert a file into a .tar archive file.
#	    The default uses the POSIX.2 pax command, which is not
#	    available on traditional UNIX systems. These typically use
#	    a tar command such as 'tar -cvhf' instead (if symbolic
#	    links are supported, then an option to follow these must be
#	    included in this command).
#	 TAR_OUT (pax -r -f)
#	    Command to use to extract a file from a .tar archive file.
#	    The default uses the POSIX.2 pax command, which is not
#	    available on traditional UNIX systems. These typically use
#	    a tar command such as 'tar -xf' instead.
#        TAR_ADD (tar -rvhf)
#           Command to use to add a file to an existing .tar archive file. The
#           default is a traditional UNIX tar command, since hme doesn't know
#           the proper POSIX.2 command options.

#  Copyright:
#     Copyright (C) 1993 Science & Engineering Research Council

#  Authors:
#     RFWS: R.F. Warren-Smith (STARLINK, RAL)
#     PDRAPER: P.W. Draper (Starlink, Durham University)
#     HME: Horst Meyerdierks (UoE, Starlink)
#     {enter_new_authors_here}

#  History:
#     12-FEB-1993 (RFWS):
#     	 Original version.
#     24-JUN-1993 (PDRAPER):
#        Adapted to CCDPACK.
#     06-JUL-1993 (PDRAPER):
#        Added Solaris and OSF1.
#     13-AUG-1993 (PDRAPER):
#        Added the ALINK environment variable.
#     29-JUN-1994 (HME):
#        Added the TAR_ADD environment variable.
#     {enter_changes_here}

#  Bugs:
#     {note_any_bugs_here}

#.

#  Export "local" definitions to the environment for use by make.
      export ALINK
      export AR_IN
      export BLD_SHR
      export CC
      export CFLAGS
      export FC
      export FFLAGS
      export LINK
      export RANLIB
      export TAR_IN
      export TAR_OUT
      export TAR_ADD

#  Check that the SYSTEM environment variable is defined.
      if test "$SYSTEM" = ""; then
         echo "mk: Please define the environment variable SYSTEM to identify"
         echo "    your computer system (the prologue in the mk script file"
         echo "    contains more information if you require it)."

#  If OK, test for each recognised system.
      else
         case "$SYSTEM" in

#  DEC Alpha:
#  =========
#  DEC Alpha machines running OSF1.
#  -------------------------------
            alpha_OSF1)
               FC='f77'
               FFLAGS='-O -w'
               CC='cc'
               LINK='ln -s'
               RANLIB='ranlib'
               TAR_IN='tar -cvhf'
               TAR_OUT='tar -xmf'
               TAR_ADD='tar -rvhf'
               echo "mk: Environment variables defined for $SYSTEM system"
	       ;;

#  DECstations:
#  ===========
#  DECstations running Ultrix.
            mips)
               FC='f77'
               FFLAGS='-O -w'
               CC='c89'
               LINK='ln -s'
               RANLIB='ranlib'
               TAR_IN='tar -cvhf'
               TAR_OUT='tar -xmf'
               TAR_ADD='tar -rvhf'
               echo "mk: Environment variables defined for $SYSTEM system"
	       ;;

#  SUN4 systems:
#  ============
#  SUN Sparcstations running SunOS 4.x.
#  -----------------------------------
            sun4)
               FC='f77'
               FFLAGS='-O -w -temp=.'
               CC='gcc'
               LINK='ln -s'
               RANLIB='ranlib'
               TAR_IN='tar -cvhf'
               TAR_OUT='tar -xmf'
               TAR_ADD='tar -rvhf'
               echo "mk: Environment variables defined for $SYSTEM system"
               ;;

#  SUN Sparcstations running SunOS 5.x (Solaris).
#  ---------------------------------------------
            sun4_Solaris)
               FC='f77'
               FFLAGS='-g -w -temp=.'
               CC='cc'
               LINK='ln -s'
               TAR_IN='tar -cvhf'
               TAR_OUT='tar -xmf'
               TAR_ADD='tar -rvhf'
               echo "mk: Environment variables defined for $SYSTEM system"
               ;;

#  Issue a warning if SYSTEM is not recognised.
	    *)
               echo "mk: WARNING: value of SYSTEM = $SYSTEM not recognised..."
               echo "             ...assuming default system characteristics"
               ;;
         esac

#  Invoke make with the appropriate environment variables set to over-ride
#  default macros defined in the makefile.
         echo make -e $*
         make -e $*
      fi

#  End of script.
