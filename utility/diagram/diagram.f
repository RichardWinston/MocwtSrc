c EXAMPLE OF OUTPUT (looks better if you choose IBM PC line graphics):

c      +---------------- subroutine a(x)             |   1
c      |+--------------- do i=1,5                    |   2
c      ||+---------------- if(i/2*2.eq.i)then        |   3
c      |||                   x=x*i                   |   4
c      ||+---------------- else                      |   5
c      |||                   x=x/i                   |   6
c      ||+---------------- endif                     |   7
c      |+--------------- enddo                       |   8
c      +---------------- end                         |   9

c Diagrams FORTRAN if-else-elseif-endif, do-enddo and case constructs,
c  start and end of routines, type definitions, modules and interfaces;
c  puts a * next to goto, return, cycle, exit, stop, end= and err=.

c Program by Mitchell R Grunes, ATSC/NRL (grunes@nrlvax.nrl.navy.mil).
c Revision date: 12/1/95.
c If you find it useful, or find a problem, please send me e-mail.

c This program was written in FORTRAN, the One True Language.

c This was written in Fortran 77 (with common extensions) for
c  portability.  It should also compile under Fortran 90 and Fortran 95,
c  provided you tell the compiler it is in card format.

c It can be confused if an INCLUDE block contains a structure that
c  begins inside and ends outside (or vice-versa).

c I hope this works for you, but bear in mind that nothing short of
c  a full-fledged language parser could really do the job.  Perhaps
c  worth about what you paid for it.    (-:

c Versions: To diagram Fortran:     diagramf.for
c                      IDL/PV-WAVE: diagrami.for
c                      C:           diagramc.for
c MS-DOS procedures to call above programs without asking so many
c  questions, append output to file diagram.out:
c                      Fortran:     diagramf.bat (card format)
c                                   diagram9.bat (free format)
c                      IDL/PV-WAVE: diagrami.bat
c                      C:           diagramc.bat
c Similar Unix csh procedures:
c                      Fortran:     diagramf.sh  (card format)
c                                   diagram9.sh  (free format)
c                      IDL/PV-WAVE: diagrami.sh
c                      C:           diagramc.sh

        program diagramf                        ! Diagrammer for Fortran
        character*80 filnam,filnam2

        print*,'FORTRAN source filename?'
        read(*,'(a80)')filnam
        print*,filnam

        print*,'Output file (blank=screen)?'
        read(*,'(a80)')filnam2
        print*,filnam2

        print*,'Column in which to write line #''s ',
     &   '(0 for none; 67 for 80 col screen; 73 to show card format):'
        read*,LCol
        print*,LCol

        print*,'Embed include files (0=no):'
        read*,iembed
        print*,iembed
        print*,' '
        print*,'0=Card format (cols 1-6 special, warnings past 72)'
        print*,'1=Free format'
        print*,'2=Card format (same as 0, ignore cols past 72)'
        print*,'Format #:'
        read*,ifree
        print*,ifree

        print*,'Use IBM PC graphics characters (0=no):'
        read*,igraphics
        print*,igraphics

        call diagram(filnam,filnam2,LCol,iembed,ifree,igraphics)
        end
c-----------------------------------------------------------------------
        subroutine diagram(filnam,filnam2,LCol,iembed,ifree,igraphics)
c Program by Mitchell R Grunes, ATSC/NRL (grunes@nrlvax.nrl.navy.mil).
        character*80 filnam,filnam2
        character*160 a,b,AfterSemi
        character*5 form
        character*8 fm
        character*1 c,c2
        logical find
        external find
        common iCol,iCol1
        character*10 label(100)
        logical fout

c Symbols which will mark block actions:
        character*1 BlockBegin    (2) /'+','+'/  ! Start of block
        character*1 BlockEnd      (2) /'+','+'/  ! End of block
        character*1 BlockElse     (2) /'+','+'/  ! Else construct
        character*1 BlockContinue (2) /'|','|'/  ! Block continues w/o change
        character*1 BlockHoriz    (2) /'-','-'/  ! Horizontal to start of line
c Same, but allows horizontal line to continue through:
        character*1 BlockBeginH   (2) /'+','+'/  ! Start of block
        character*1 BlockEndH     (2) /'+','+'/  ! End of block
        character*1 BlockElseH    (2) /'+','+'/  ! Else construct

        if(iGraphics.ne.0)then
          iGraphics=1

          BlockBegin   (1)=char(218)            ! (1)=normal
          BlockEnd     (1)=char(192)
          BlockElse    (1)=char(195)
          BlockContinue(1)=char(179)
          BlockHoriz   (1)=char(196)
          BlockBeginH  (1)=char(194)
          BlockEndH    (1)=char(193)
          BlockElseH   (1)=char(197)

          BlockBegin   (2)=char(214)            ! (2)=DO/FOR loops (doubled)
          BlockEnd     (2)=char(211)            ! (not yet used)
          BlockEnd     (2)=char(211)
          BlockElse    (2)=char(199)
          BlockContinue(2)=char(186)
          BlockHoriz   (2)=char(196)
          BlockBeginH  (2)=char(209)
          BlockEndH    (2)=char(208)
          BlockElseH   (2)=char(215)
        endif

        open(1,file=filnam,status='old')
        fout=filnam2.gt.' '
        if(fout)open(2,file=filnam2,status='unknown')
                                                ! ASCII 12 is a form feed
        if(fout)write(2,*)char(12),
     &   '=============--',filnam(1:LenA(filnam)),'--============='

        if(fout)     write(2,'(11x,a50,a49,/)')  ! Write column header
     &   '....,....1....,....2....,....3....,....4....,....5',
     &   '....,....6....,....7....,....8....,....9....,....'
        if(.not.fout)write(*,'(11x,a50,a49,/)')' ',
     &   '....,....1....,....2....,....3....,....4....,....5',
     &   '....,....6....,....7....,....8....,....9....,....'

        i1=0                                    ! # of nest levels before
                                                !  current line
        i2=0                                    ! # of nest levels on
                                                !  current line
        i3=0                                    ! # of nest levels after
                                                !  current line
        i4=0                                    ! not 0 to flag start or end
                                                !  of block
        InSub=0                                 ! Inside a subroutine,
                                                !  function or mainline
        InMod=0                                 ! Inside module or
                                                !   contains
        nMain=0                                 ! no mainline program yet
        InElse=0                                ! Found elseif, but not then
        nlabel=0                                ! # of labels for do loop
                                                !  ends
        iAlphaNum=0                             ! Last char of line is
                                                !  alpha-numeric
        iOldContinue=0                          ! next line not continued line
        nline=0
        iunit=1
10      a=' '
        read(iunit,'(a160)',end=99)a
        nline=nline+1
        fm=' '
        write(fm,'(i5)')nline
        form=fm

        if(a(1:1).eq.char(12))then
          if(fout)write(2,'(a1,:)')char(12)
          if(.not.fout)print*,'------------FORM FEED------------'
          b=a(2:160)
          a=b
        endif

        b=' '                                   ! Turn tabs to spaces
        j=1
        do i=1,LenA(a)
          if(a(i:i).eq.char(9))then
            j=(j-1)/8*8+8+1
          elseif(j.le.160)then
            b(j:j)=a(i:i)
            j=j+1
          endif
        enddo

        a=' '                                   ! Pre-processed output
        i=1                                     ! Basic pre-processing
        j=1
        i72flag=0                               ! nothing over column 72
                                                ! yet
        iOldAlphaNum=iAlphaNum                  ! last line ended in
                                                ! alpha-numeric?
        iAlphaNum=0
        iContinue=iOldContinue                  ! This line continued line?
        if(find(b,'&',2,0))iContinue=1          !  will be changed to 2 after
                                                !  first non/blank.
        if(i.eq.6.and.c.ne.' '.and.ifree.ne.1)iContinue=1
        if(iContinue.eq.0)then
          iquote=0                              ! no ' yet
          idquote=0                             ! no " yet
        endif
        j=1
                                                ! comment line
        if((b(1:1).eq.'c'.or.b(1:1).eq.'C').and.ifree.ne.1)goto 15


        do i=1,LenA(b)
          c=b(i:i)
                                                ! handle upper case
          if(c.ge.'A'.and.c.le.'Z')c=char(ichar(c)+32)
                                                ! ASCII 33 is '!'
          if(c.eq.char(33).and.iquote.eq.0.and.idquote.eq.0)goto 15

          if(i.gt.72.and.c.ne.' ')then
            if(ifree.eq.0.and.i72flag.eq.0)then
              i72flag=1
              PRINT*,'***WARNING--PAST COLUMN 72 at line',form
              if(fout)print*,b
              print*,char(7)
            elseif(ifree.eq.2)then
              c=' '
            endif
          endif

          if(c.eq.''''.and.(i.ne.6.or.ifree.ne.0).and.idquote.eq.0)
     &     iquote=1-iquote
          if(c.eq.'"' .and.(i.ne.6.or.ifree.ne.0).and.iquote .eq.0)
     &     idquote=1-idquote
          if(iquote.eq.1)then
            if(find(a,'include ',2,0).and.iembed.ne.0)then
              iquote=0
              idquote=0
            endif
          endif
          if(iquote.ne.0.or.idquote.ne.0)c=' '
          if(j.gt.1)then                        ! (kill multiple spaces,
                                                !  and spaces around =)
            c2=a(j-1:j-1)
            if(c.eq.' '.and.c2.eq.' ')j=j-1
            if(c.eq.'='.and.c2.eq.' ')j=j-1
            if(c.eq.' '.and.c2.eq.'=')j=j-1
            if(c.eq.' '.and.c2.eq.'=')c='='
          endif
                                                ! Look for
                                                ! identifiers that wrap
                                                ! around lines.
          if((i.gt.6.or.ifree.ne.0).and.c.ne.' '.and.c.ne.'&')then
            iAlphaNum=0
            if((c.ge.'a'.and.c.le.'z').or.
     &       (c.ge.'0'.and.c.le.'9'))then
              iAlphaNum=1
              if(iContinue.eq.1)then
                if(iOldAlphaNum.ne.0)then
                  PRINT*,'***POSSIBLE SPLIT IDENTIFIER across line',form
                  print*,char(7)
                endif
              endif
            endif
            iContinue=2
          endif

          if(j.le.160) a(j:j)=c
          j=j+1
        enddo

15      iOldContinue=0
        if(a(LenA(a):LenA(a)).eq.'&')iOldContinue=1

        i2=i1
        i3=i1
        i4=0
        igoto=0                                 ! no goto on line
        Main1=0                                 ! (Not mainline)
                                                ! Possible mainline start

16      AfterSemi=' '                           ! Break line at semicolons
        if(find(a,';',0,160-1))then
          AfterSemi='      '//a(icol:160)
          a=a(1:icol1-1)
        endif

        if(a.ne.' '.and.InSub.eq.0.and.InMod.eq.0)Main1=1
                                                ! Mark various types of jump
        if(find(a,'go to',8+64,0).or.find(a,'goto',8+64,0).or.
     &     find(a,'end=',16,0)   .or.find(a,'err=',16,0)  .or.
     &     find(a,'return',8+64,0).or.find(a,'cycle ',8,0).or.
     &     find(a,'exit ',8,0)   .or.find(a,'stop ',8,0))
     &   igoto=1

        if(find(a,')1',64,0).or.find(a,')2',64,0).or.
     &     find(a,')3',64,0).or.find(a,')4',64,0).or.
     &     find(a,')5',64,0).or.find(a,')6',64,0).or.
     &     find(a,')7',64,0).or.find(a,')8',64,0).or.
     &     find(a,')9',64,0))
     &   igoto=1

        if(find(a,') 1',64,0).or.find(a,') 2',64,0).or.
     &     find(a,') 3',64,0).or.find(a,') 4',64,0).or.
     &     find(a,') 5',64,0).or.find(a,') 6',64,0).or.
     &     find(a,') 7',64,0).or.find(a,') 8',64,0).or.
     &     find(a,') 9',64,0))
     &   igoto=1

        if(find(a,'::',0,0))then              ! To distinguish
          iDeclare=iCol                         !  declarations from
                                                !  keywords
        else
          iDeclare=999
        endif

        if(find(a,'include ''',2,0).and.iembed.ne.0)then
          filnam=a(iCol:160)
          if(.not.find(filnam,'''',0,0))goto 20
          filnam(iCol-1:80)=' '
          if(fout)print*,'including file ',filnam(1:50)
          close(3)
          open(3,file=filnam,status='old',err=17)
          iunit=3
          nlinesave=nline
          nline=0
          i2=i2+1
          i3=i3+1
          goto 20
17        PRINT*,'***WARNING--Missing include file***'
          print*,char(7)
        elseif(find(a,'end module ',2,0).or.
     &         find(a,'endmodule ',2,0).or.
     &         find(a,'end interface',2,0).or.
     &         find(a,'endinterface',2,0).or.
     &         find(a,'end type ',2,0).or.
     &         find(a,'endtype ',2,0))then
          i3=i3-1
          InMod=InMod-1
          if(find(a,'endmodule ',2,0).or.
     &       find(a,'end module ',2,0))then
            InMod=0
            if(InSub.gt.0.or.i3.ne.0)then
              PRINT*,'***ERROR--INVALID DIAGRAMMING INDEX line',form
              if(fout)WRITE(2,*)
     &         '***ERROR--INVALID DIAGRAMMING INDEX!***'
              if(fout)print*,b
              print*,char(7)
            endif
          endif
          InElse=0
        elseif(find(a,'enddo ',256,0).or.
     &         find(a,'end do  ',256,0))then
          i3=i3-1
          nlabel=max(0,nlabel-1)
          InElse=0
        elseif(find(a,'endif  ',256,0).or.
     &         find(a,'end if  ',256,0).or.
     &         find(a,'endselect ',256,0).or.
     &         find(a,'end select ',256,0).or.
     &         find(a,'endforall ',256,0).or.
     &         find(a,'end forall ',256,0).or.
     &         find(a,'endforall ',256,0).or.
     &         find(a,'end where ',256,0))then
          i3=i3-1
          InElse=0
        elseif(find(a,'end  ',256,0).or.
     &         find(a,'end function ',256,0).or.
     &         find(a,'endfunction ',256,0).or.
     &         find(a,'end subroutine ',256,0).or.
     &         find(a,'endsubroutine ',256,0).or.
     &         find(a,'end program ',256,0).or.
     &         find(a,'endprogram ',256,0).or.
     &         find(a,'end block',256,0).or.
     &         find(a,'endblock',256,0))then
          i3=i3-1
          InSub=InSub-1
          if(InSub.lt.0.or.(InSub.gt.0.and.InMod.le.0))then
            if(InSub.lt.0.and.InMod.gt.0.and.find(a,'end  ',256,0))then
              InSub=0
              InMod=InMod-1
            else
              PRINT*,'***ERROR--INVALID DIAGRAMMING INDEX line',form
              if(fout)
     &         WRITE(2,*)'***ERROR--INVALID DIAGRAMMING INDEX!***'
              if(fout)print*,b
              print*,char(7)
            endif
          endif
          if(i3.eq.0)InSub=0
          InElse=0
        elseif(find(a,'elseif',128+256,0).or.
     &         find(a,'else if',128+256,0))then
          i4=max(i4,1)
          InElse=0
          if(.not.find(a,'then  ',8,0))InElse=1
        elseif(find(a,'then  ',8,0))then
          i2=i2+1
          if(InElse.eq.0)i3=i3+1
          InElse=0
        elseif( find(a,'selectcase',256,0).or.
     &          find(a,'select case',256,0))then
          i2=i2+1
          i3=i3+1
          i4=max(i4,1)
          InElse=0
        elseif(find(a,'else  ',256,0).or.
     &         find(a,'entry ',4,0).or.
     &         find(a,'case ',256,0).or.
     &         find(a,'case(',256,0).or.
     &         find(a,'contains ',2,0).or.
     &         find(a,'elsewhere  ',256,0).or.
     &         find(a,'else where  ',256,0))then
          i4=max(i4,1)
          InElse=0
          if(find(a,'contains ',2,0))then
            if(fout)print*,'Line ',form,' ',b(1:LenA(b))
            InMod=InMod+1
          endif
        elseif( find(a,'selectcase',256,0).or.
     &          find(a,'select case',256,0).or.
     &          find(a,'for all (',256,0).or.
     &          find(a,'forall (',256,0).or.
     &          find(a,'for all(',256,0).or.
     &          find(a,'forall(',256,0).or.
     &          find(a,'where (',256,0).or.
     &          find(a,'where(',256,0))then
          i2=i2+1
          i3=i3+1
          InElse=0
        elseif((find(a,'module ',2,iDeclare).and.
     &        .not.find(a,'module procedure',2,iDeclare)).or.
     &         find(a,'interface ',2,iDeclare).or.
     &        (find(a,'type ',2,iDeclare).and.
     &        .not.find(a,'(',0,iDeclare)))then
          if(fout)print*,'Line ',form,' ',b(1:LenA(b))
          i2=i2+1
          i3=i3+1
          Main1=0
          if(find(a,'module ',2,iDeclare).and.InMod.ne.0)then
            PRINT*,'***ERROR--NESTED MODULES***'
            if(fout)WRITE(2,*)'***NESTED MODULES***'
            if(fout)print*,b
            print*,char(7)
          endif
          InMod=InMod+1
          InElse=0
        elseif(find(a,'do while',128+256,0).or.
     &         find(a,'dowhile',128+256,0))then
          i2=i2+1
          i3=i3+1
          nlabel=min(100,nlabel+1)
          label(nlabel)='####'
          InElse=0
        elseif(find(a,' do ',256,0).or.
     &   (ifree.ne.0.and.a(1:3).eq.'do '))then
         if(ifree.ne.0.and.a(1:3).eq.'do ')iCol=4
          if(iCol1.lt.7.or.a(7:max(7,iCol1)).eq.' '.or.
     &     (ifree.ne.0.and.a(1:3).eq.'do '))then
            i2=i2+1
            i3=i3+1
            iCol2=iCol
            dowhile(iCol2.lt.160.and.a(iCol2:iCol2).ge.'0'.and.
     &       a(iCol2:iCol2).le.'9')
              iCol2=iCol2+1
            enddo
            iCol2=iCol2-1
            nlabel=min(100,nlabel+1)
            if(iCol2.ge.iCol)then
              label(nlabel)=a(iCol:iCol2)
            else
              label(nlabel)='####'
            endif
          endif
          InElse=0
        elseif(find(a,': do ',0,0).or.find(a,':do ',0,0))then
          i2=i2+1
          i3=i3+1
          InElse=0
        elseif(find(a,'function ',4,iDeclare).or.
     &         find(a,'subroutine ',4,iDeclare).or.
     &         find(a,'program ',2,iDeclare) .or.
     &         find(a,'block data ',2,iDeclare).or.
     &         find(a,'blockdata ',2,iDeclare))then
          if(fout)print*,'Line ',form,' ',b(1:LenA(b))
          if(InSub.ne.0.and.InMod.eq.0)then
            PRINT*,'***ERROR--ROUTINE INSIDE ROUTINE***'
            if(fout)WRITE(2,*)'***ERROR--ROUTINE INSIDE ROUTINE***'
            if(fout)print*,b
            print*,char(7)
          endif
          Main1=0
          InSub=InSub+1
          i2=i2+1
          i3=i3+1
          if(InSub.eq.1.and.i3.ne.1.and.InMod.le.0)then
            PRINT*,'***ERROR--INVALID DIAGRAMMING INDEX line',form
            if(fout)
     &       WRITE(2,*)'***ERROR--INVALID DIAGRAMMING INDEX!***'
            if(fout)print*,b
            print*,char(7)
            i3=1
          endif
          InElse=0
        endif

20      if(Main1.ne.0)then                      ! Was start of mainline
          if(fout)print*,'Line ',form,' ',b(1:LenA(b))
          if(nMain.gt.0)then
            PRINT*,'***ERROR--TOO MANY MAINLINES***'
            if(fout)WRITE(2,*)'***ERROR--TOO MANY MAINLINES!***'
            if(fout)print*,b
            print*,char(7)
          endif
          InSub=InSub+1
          nMain=nMain+1
          i2=i2+1
          i3=i3+1
        endif

21      if(b(1:5).ne.' '.or.ifree.ne.0)then     ! Search for DO labels
          iend=1
          dowhile(iend.lt.160.and.(b(iend:iend).eq.' '.or.
     &     (b(iend:iend).ge.'0'.and.b(iend:iend).le.'9')))
            iend=iend+1
          enddo
          iend=iend-1
          if(iend.ge.1.and.b(1:max(1,iend)).ne.' ')then
            do i=1,nlabel
              j=nlabel+1-i                      !  (in reverse order)
              if(find(b(1:iend),label(j)(1:LenA(label(j))),1,0))then
                i3=i3-1
                nlabel=max(0,j-1)
                goto 21
              endif
            enddo
          endif
        endif

        if(AfterSemi.ne.' ')then
          a=AfterSemi
          goto 16
        endif

        a=' '
        if(i1.lt.0.or.i2.lt.0.or.i3.lt.0.or.i4.lt.0)then
          PRINT*,'***ERROR--INVALID DIAGRAMMING INDEX line',form
          if(fout)WRITE(2,*)'***ERROR--INVALID DIAGRAMMING INDEX!***'
          if(fout)print*,b
          print*,char(7)
          i1=max(i1,0)
          i2=max(i2,0)
          i3=max(i3,0)
          i4=max(i4,0)
        endif

        i2=max(i1,i3)                           ! # of nests on current line
        i4=max(i4,iabs(i3-i1))                  ! not 0, to flag start or
                                                !  end of block

        iBlock=1                                ! For the present version.

        a=' '                                   ! Leave space for diagram
        a(12:160)=b                             !  (must match column header)

        LastUse=1                               ! Last usable diagram col
        dowhile(LastUse.lt.160.and.a(LastUse:LastUse).eq.' ')
          LastUse=LastUse+1
        enddo
        LastUse=LastUse-2

        if(igoto.ne.0)a(1:1)='*'                ! Place * next to jumps

        if(i2.gt.0)then                         ! Draw one vertical line per
          do i=2,min(i2+1,LastUse)              !  nest level.
            a(i:i)=BlockContinue(iBlock)
          enddo
        endif

        if(i4.ne.0)then                         ! Draw horizontal lines inward
          do i=i2+2,LastUse                     !  from above.
            a(i:i)=BlockHoriz(iBlock)
          enddo
        endif

        do i=0,i4-1                             ! May need to replace some
                                                !  vertical lines with
          c=              BlockElse(iBlock)     !       else  symbol
          if(i1+i.lt.i3)c=BlockBegin(iBlock)    !  or   begin symbol
          if(i1+i.gt.i3)c=BlockEnd   (iBlock)   !  or   end   symbol
          j=max(2,min(LastUse,i2+1-i))
          a(j:j)=c
          if(a(j+1:j+1).eq.BlockElse  (iBlock)) ! Continue horizontal lines
     &       a(j+1:j+1)  = BlockElseH (iBlock)
          if(a(j+1:j+1).eq.BlockBegin (iBlock))
     &       a(j+1:j+1)  = BlockBeginH(iBlock)
          if(a(j+1:j+1).eq.BlockEnd   (iBlock))
     &       a(j+1:j+1)  = BlockEndH  (iBlock)
        enddo

        if(LCol.gt.0.and.a(max(1,LCol+11):160).eq.' ')then       ! line #
          if(form(1:1).eq.' ')form(1:1)=BlockContinue(iBlock)
          a(LCol+11:160)=form
        endif

        n=LenA(a)                               ! Output diagrammed line
        if(fout)     write(2,'(80a1,80a1)')(a(i:i),i=1,n)
        if(.not.fout)write(*,'(1x,80a1,80a1)')(a(i:i),i=1,n)

        i1=i3
        goto 10
99      if(iunit.eq.3)then
          iunit=1
          i1=i1-1
          close(3)
          nline=nlinesave
          goto 10
        endif
        if(i3.gt.0.or.InSub.ne.0)then
          PRINT*,'***WARNING--SOME NEST LEVELS LEFT HANGING AT END***'
          print*,char(7)
        endif
        end
c-----------------------------------------------------------------------
        logical function find(a,b,icond,jcol)   ! find b in a, subject
                                                !  to conditions:
                                                ! Colunn is prior to jcol
                                                !  (if jcol.ne.0)
                                                ! icond=sum of the
                                                !  following:
                                                ! 1:  Prior, if exists, must
                                                !     be blank
                                                ! 2:  Must be first non-blank
                                                ! 4:  Prior character, if
                                                !  present, must not be
                                                !  alphanumeric.
                                                ! 8:  Prior character, if
                                                !  present, must be blank
                                                !  or )
                                                ! 16: Prior character, if
                                                !  present, must be blank
                                                !  or ,
                                                ! 32: Next character not
                                                !  alphanumeric
                                                ! 64: Next character not
                                                !  alphabetic
                                                ! 128:Next character must
                                                !  be blank or (
                                                ! 256:1st non-blank,
                                                !  possibly except for
                                                !  numeric labels
                                                ! 512  Prior character, if present,
                                                !      must be blank or ) or }
                                                !      or { or ;
c Program by Mitchell R Grunes, ATSC/NRL (grunes@nrlvax.nrl.navy.mil).
c Revision date: 11/30/95.
        character*(*) a,b
        character*1   c,cNext,c2
        common iCol,iCol1
        logical result

        ii=len(a)
        jj=len(b)
        result=.false.
        jjcol=999
        if(jcol.gt.0)jjcol=jcol
        do i=1,min(ii-jj+1,jjcol)
          if(a(i:i+jj-1).eq.b)then              ! Found--Now do tests
            iCol1=i                             ! iCol1=column of item
                                                !  found
            iCol =i+jj                          ! iCol =colomn after
                                                !  item found

            c=' '
            cNext=' '
            if(iCol1.gt.1)c=a(iCol1-1:iCol1-1)
            if(iCol .le.ii)cNext=a(iCol:iCol)

            result=.true.
            if(result.and.iand(icond,1).ne.0.and.icol1.gt.1)then
              result=c.eq.' '
            endif

            if(result.and.iand(icond,2).ne.0.and.iCol1.gt.1)then
              result=a(1:iCol1-1).eq.' '
            endif

            if(result.and.iand(icond,4).ne.0)
     &       result=(c.lt.'0'.or.c.gt.'9').and.(c.lt.'a'.or.c.gt.'z')
            if(result.and.iand(icond,8).ne.0)result=c.eq.' '.or.c.eq.')'

            if(result.and.iand(icond,16).ne.0)
     &       result=c.eq.' '.or.c.eq.','

            if(result.and.iand(icond,32).ne.0)
     &       result=(cNext.lt.'0'.or.cNext.gt.'9').and.
     &              (cNext.lt.'a'.or.cNext.gt.'z')

            if(result.and.iand(icond,64).ne.0)
     &       result=(cNext.lt.'a'.or.cNext.gt.'z')

            if(result.and.iand(icond,128).ne.0)
     &       result=cNext.eq.' '.or.cNext.eq.'('

            if(result.and.iand(icond,256).ne.0.and.iCol1.gt.1)then
              do iii=1,iCol1-1
                c2=a(iii:iii)
                if((c2.lt.'0'.or.c2.gt.'9').and.c2.ne.' ')result=.false.
              enddo
            endif

            if(result.and.iand(icond,512).ne.0)result=c.eq.' '
     &       .or.c.eq.';'.or.c.eq.')'.or.c.eq.'{'.or.c.eq.'}'

            find=result
            if(result)return
          endif
        enddo
        find=result
        end
c-----------------------------------------------------------------------
        function LenA(a)                        ! Length of string, at
                                                !  least 1
c Program by Mitchell R Grunes, ATSC/NRL (grunes@nrlvax.nrl.navy.mil).
c Revision date: 11/30/95.
        character*(*) a
        n=len(a)
        dowhile(n.gt.1.and.a(n:n).eq.' ')
          n=n-1
        enddo
        LenA=n
        end
