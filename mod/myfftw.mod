	  f  Ú   k820309              12.0        *gàX                                                                                                           
       myfftw.f90 MYFFTW       	       FFTW_ESTIMATE DLC IU PI TWOPI gen@SAVEARRAY+MYUTILS gen@MAPCUT+ARRAY gen@SPIN_WEIGHT+ANAFLAT gen@GAUSSIAN_ALM+ANAFLAT                                                     
       DLC IU PI TWOPI                      @                              
       SAVEARRAY                      @                              
       MAPCUT                      @                              
       SPIN_WEIGHT GAUSSIAN_ALM                                                       u #DFT_1DARRAY    #DFT_2DARRAY                                                          u #DFT_POL_1DARRAY    #DFT_POL_2DARRAY                                                           u #DFT_ALL_1DARRAY 	   #DFT_ALL_2DARRAY 
                                                          u #DERIVEMAP_ALL    #DERIVEMAP_NTH_1D    #DERIVEMAP_NTH_2D    #DERIVEMAP_ALM_NTH_1D    #         @                                                   #SAVEARRAY_2D_DBLE%TRIM    #SAVEARRAY_2D_DBLE%SIZE    #F    #ARR                                                    TRIM                                                 SIZE           
                                                     1           
                                                    
 _             &                   &                                           #         @                                                   #SAVEARRAY_2D_CMPLX%TRIM    #SAVEARRAY_2D_CMPLX%SIZE    #F    #ARR                                                    TRIM                                                 SIZE           
                                                     1           
                                                    `             &                   &                                           #         @                                                   #SAVEARRAY_1D_DBLE%RESHAPE    #F    #ARR    #NN                                                    RESHAPE           
                                                     1           
                                                     
 \             &                                                     
                                                      [   p          & p        p            p                          #         @                                                   #SAVEARRAY_1D_CMPLX%RESHAPE    #F     #ARR !   #NN "                                                   RESHAPE           
                                                      1           
                                 !                    ^             &                                                     
                                  "                    ]   p          & p        p            p                          #         @                                #                  #MAPCUT_DBLE%RESHAPE $   #MAPCUT_DBLE%TRANSPOSE %   #IMAP &   #MM '   #OMAP (   #NN )                                              $     RESHAPE                                            %     TRANSPOSE           
                                  &                   
              &                                                     
                                  '                       p          & p        p            p                                                                     (                   
               &                                                     
                                  )                       p          & p        p            p                          #         @                                *                  #MAPCUT_CMPLX%DBLE +   #MAPCUT_CMPLX%AIMAG ,   #IMAP -   #MM .   #OMAP /   #NN 0                                              +     DBLE                                            ,     AIMAG           
                                 -                                 &                                                     
                                  .                       p          & p        p            p                                                                    /                                  &                                                     
                                  0                       p          & p        p            p                          #         @                                1                  #SPIN_WEIGHT_1DARRAY%ABS 2   #SPIN_WEIGHT_1DARRAY%DBLE 3   #SW 4   #NN 5   #D 6   #TRANS 7   #S 8                                              2     ABS                                            3     DBLE                                           4                                  &                                                     
                                  5                       p          p            p                                    
                                  6                   
    p          p            p                                    
                                  7                     
                                  8           #         @                                9                  #SPIN_WEIGHT_2DARRAY%ABS :   #SPIN_WEIGHT_2DARRAY%DBLE ;   #SW <   #NN =   #D >   #TRANS ?   #S @                                              :     ABS                                            ;     DBLE                                           <                                  &                   &                                                     
                                  =                       p          p            p                                    
                                  >                   
    p          p            p                                    
                                  ?                     
                                  @           #         @                                A              	    #GAUSSIAN_ALM_1DARRAY%PRESENT B   #GAUSSIAN_ALM_1DARRAY%ABS C   #GAUSSIAN_ALM_1DARRAY%INT D   #GAUSSIAN_ALM_1DARRAY%DBLE E   #GAUSSIAN_ALM_1DARRAY%DSQRT F   #GAUSSIAN_ALM_1DARRAY%SIZE G   #GAUSSIAN_ALM_1DARRAY%MOD H   #GAUSSIAN_ALM_1DARRAY%CMPLX I   #GAUSSIAN_ALM_1DARRAY%CONJG J   #NN K   #D L   #IL M   #ALM N   #CL O   #FIX P                                              B     PRESENT                                            C     ABS                                            D     INT                                            E     DBLE                                            F     DSQRT                                            G     SIZE                                            H     MOD                                            I     CMPLX                                            J     CONJG           
                                  K                    l   p          p            p                                    
                                  L                   
 n   p          p            p                                    
                                 M                    k   p          p            p                                    
                                N                    o              &                                                     
                                 O                   
 m             &                                                     
                                 P           #         @                                 Q              	    #GAUSSIAN_ALM_2DARRAY%PRESENT R   #GAUSSIAN_ALM_2DARRAY%ABS S   #GAUSSIAN_ALM_2DARRAY%INT T   #GAUSSIAN_ALM_2DARRAY%DBLE U   #GAUSSIAN_ALM_2DARRAY%DSQRT V   #GAUSSIAN_ALM_2DARRAY%SIZE W   #GAUSSIAN_ALM_2DARRAY%MOD X   #GAUSSIAN_ALM_2DARRAY%CMPLX Y   #GAUSSIAN_ALM_2DARRAY%CONJG Z   #NN [   #D \   #IL ]   #ALM ^   #CL _   #FIX `                                              R     PRESENT                                            S     ABS                                            T     INT                                            U     DBLE                                            V     DSQRT                                            W     SIZE                                            X     MOD                                            Y     CMPLX                                            Z     CONJG           
                                  [                    c   p          p            p                                    
                                  \                   
 e   p          p            p                                    
                                 ]                    b   p          p            p                                    
                                ^                    f              &                   &                                                     
                                 _                   
 d             &                                                     
                                 `           #         @      X                                             #DFT_1DARRAY%ABS a   #DFT_1DARRAY%SUM b   #DFT_1DARRAY%DBLE c   #MAP d   #NN e   #D f   #TRANS g                                              a     ABS                                            b     SUM                                            c     DBLE           
D @                              d                                  &                                                     
@ @                               e                    
   p          & p        p            p                                    
                                  f                   
    p          & p        p            p                                    
@ @                               g           #         @      X                                             #DFT_2DARRAY%ABS h   #DFT_2DARRAY%SUM i   #DFT_2DARRAY%DBLE j   #MAP k   #NN l   #D m   #TRANS n                                              h     ABS                                            i     SUM                                            j     DBLE           
D @                              k                                  &                   &                                                    
@ @                               l                       p          p            p                                    
                                  m                   
    p          p            p                                    
@ @                               n           #         @      X                                             #DFT_POL_1DARRAY%CONJG o   #DFT_POL_1DARRAY%AIMAG p   #DFT_POL_1DARRAY%DBLE q   #QU r   #NN s   #D t   #EB u   #TRANS v                                              o     CONJG                                            p     AIMAG                                            q     DBLE           
D                                 r                   
               &                   &                                                     
@ @                               s                       p          p            p                                    
@ @                               t                   
    p          p            p                                    
D                                u                                  &                   &                                                     
                                  v           #         @      X                                              #DFT_POL_2DARRAY%CONJG w   #DFT_POL_2DARRAY%AIMAG x   #DFT_POL_2DARRAY%DBLE y   #QU z   #NN {   #D |   #EB }   #TRANS ~                                              w     CONJG                                            x     AIMAG                                            y     DBLE           
D                                 z                   
               &                   &                   &                                                     
@ @                               {                       p          p            p                                    
@ @                               |                   
    p          p            p                                    
D                                }                                  &                   &                   &                                                     
                                  ~           #         @      X                            	                  #DFT_ALL_1DARRAY%CONJG    #DFT_ALL_1DARRAY%AIMAG    #DFT_ALL_1DARRAY%DBLE    #TQU    #NN    #D    #TEB    #TRANS                                                    CONJG                                                 AIMAG                                                 DBLE           
D                                                    
               &                   &                                                     
@ @                                                      p          p            p                                    
@ @                                                  
    p          p            p                                    
D                                                                  &                   &                                                     
                                             #         @      X                            
                  #DFT_ALL_2DARRAY%CONJG    #DFT_ALL_2DARRAY%AIMAG    #DFT_ALL_2DARRAY%DBLE    #TQU    #NN    #D    #TEB    #TRANS                                                    CONJG                                                 AIMAG                                                 DBLE           
D                                                    
 $              &                   &                   &                                                     
@ @                                                   "   p          p            p                                    
@ @                                                  
 #   p          p            p                                    
D                                                    %              &                   &                   &                                                     
                                             #         @      X                                              #DERIVEMAP_ALL%SIZE    #DERIVEMAP_ALL%DBLE    #D    #MAP    #DERIVMAP                                                    SIZE                                                 DBLE           
@ @                                                  
 i   p          p            p                                    
 @                                                  
 j             &                   &                                                     D@                                                  
 k              &                   &                   &                   &                                           #         @      X                                              #DERIVEMAP_NTH_1D%CMPLX    #DERIVEMAP_NTH_1D%DBLE    #NN    #D    #MAP    #DMAP    #NTH                                                    CMPLX                                                 DBLE           
@ @                                                   X   p          p            p                                    
@ @                                                  
 Y   p          p            p                                    
  @                                                  
 Z             &                                                     D                                                    
 [              &                   &                                                     
                                             #         @      X                                              #DERIVEMAP_NTH_2D%SIZE    #DERIVEMAP_NTH_2D%DBLE    #D    #MAP    #DERIVMAP    #NTH                                                     SIZE                                                 DBLE           
@ @                                                  
 ^   p          p            p                                    
 @                                                  
 _             &                   &                                                     D                                                    
 `              &                   &                   &                                                     
                                              #         @      X                                              #DERIVEMAP_ALM_NTH_1D%DBLE ¡   #NN ¢   #D £   #ALM ¤   #DMAP ¥   #NTH ¦                                              ¡     DBLE           
@ @                               ¢                    d   p          p            p                                    
@ @                               £                   
 e   p          p            p                                    
                                 ¤                    f             &                                                     D                                 ¥                   
 g              &                   &                                                     
                                  ¦           #         @                                  §                  #QU2EB%DBLE ¨   #NN ©   #D ª   #QU «   #EB ¬   #TRANS ­                                              ¨     DBLE           
                                  ©                       p          p            p                                    
                                  ª                   
    p          p            p                                    
D                                «                                  &                   &                                                     
D                                ¬                                  &                   &                                                     
                                  ­           #         @                                  ®                  #PUREEB%PRESENT ¯   #PUREEB%DBLE °   #QU ±   #NN ²   #D ³   #EB ´   #W µ   #WD ¶                                              ¯     PRESENT                                            °     DBLE           
                                  ±                   
 ,             &                   &                                                     
@ @                               ²                    )   p          & p        p            p                                    
@ @                               ³                   
 *   p          & p        p            p                                    D                                ´                    .              &                   &                                                     
  @                               µ                   
 +             &                                                     
 @                               ¶                   
 -             &                   &                                           #         @                                 ·                  #ARRAY_DERIV%DBLE ¸   #NN ¹   #D º   #W »   #WX ¼   #WY ½   #WXX ¾   #WXY ¿   #WYY À                                              ¸     DBLE           
@ @                               ¹                    =   p          p            p                                    
@ @                               º                   
 >   p          p            p                                    
                                  »                   
 ?             &                                                     
D @                              ¼                    @              &                                                     
D @                              ½                    A              &                                                     
D @                              ¾                    B              &                                                     
D @                              ¿                    C              &                                                     
D @                              À                    D              &                                           #         @                                  Á                  #GAUSSIAN_MAP%PRESENT Â   #NN Ã   #MM Ä   #D Å   #EL Æ   #CL Ç   #MAP È   #FIX É                                              Â     PRESENT           
@ @                               Ã                    F   p          & p        p            p                                    
@ @                               Ä                    G   p          & p        p            p                                    
@ @                               Å                   
 I   p          & p        p            p                                    
@ @                               Æ                    H   p          & p        p            p                                    
@ @                               Ç                   
 J             &                                                     D @                              È                    K              &                                                     
B @                               É           #         @                                  Ê               	   #GAUSSIAN_MAP_POL%PRESENT Ë   #NN Ì   #MM Í   #D Î   #EL Ï   #EE Ð   #BB Ñ   #QU Ò   #P Ó   #FIX Ô                                              Ë     PRESENT           
@ @                               Ì                    N   p          & p        p            p                                    
@ @                               Í                    O   p          & p        p            p                                    
@ @                               Î                   
 Q   p          & p        p            p                                    
@ @                               Ï                    P   p          & p        p            p                                    
B @                               Ð                   
 R             &                                                     
B @                               Ñ                   
 S             &                                                     F @                               Ò                   
 T              &                   &                                                     F @                              Ó                    U              &                                                     
B @                               Ô                        fn#fn    º      b   uapp(MYFFTW    @  P   J  MYCONST      J   J  MYUTILS    Ú  G   J  ARRAY    !  Y   J  ANAFLAT    z  b       gen@DFT    Ü  j       gen@DFT_POL    F  j       gen@DFT_ALL    °         gen@DERIVEMAP *   I         SAVEARRAY_2D_DBLE+MYUTILS 4   Ù  =      SAVEARRAY_2D_DBLE%TRIM+MYUTILS=TRIM 4     =      SAVEARRAY_2D_DBLE%SIZE+MYUTILS=SIZE ,   S  L   a   SAVEARRAY_2D_DBLE%F+MYUTILS .     ¤   a   SAVEARRAY_2D_DBLE%ARR+MYUTILS +   C         SAVEARRAY_2D_CMPLX+MYUTILS 5   Õ  =      SAVEARRAY_2D_CMPLX%TRIM+MYUTILS=TRIM 5     =      SAVEARRAY_2D_CMPLX%SIZE+MYUTILS=SIZE -   O  L   a   SAVEARRAY_2D_CMPLX%F+MYUTILS /     ¤   a   SAVEARRAY_2D_CMPLX%ARR+MYUTILS *   ?         SAVEARRAY_1D_DBLE+MYUTILS :   ¾  @      SAVEARRAY_1D_DBLE%RESHAPE+MYUTILS=RESHAPE ,   þ  L   a   SAVEARRAY_1D_DBLE%F+MYUTILS .   J	     a   SAVEARRAY_1D_DBLE%ARR+MYUTILS -   Ö	  ¤   a   SAVEARRAY_1D_DBLE%NN+MYUTILS +   z
         SAVEARRAY_1D_CMPLX+MYUTILS ;   ú
  @      SAVEARRAY_1D_CMPLX%RESHAPE+MYUTILS=RESHAPE -   :  L   a   SAVEARRAY_1D_CMPLX%F+MYUTILS /        a   SAVEARRAY_1D_CMPLX%ARR+MYUTILS .     ¤   a   SAVEARRAY_1D_CMPLX%NN+MYUTILS "   ¶          MAPCUT_DBLE+ARRAY 2   V  @      MAPCUT_DBLE%RESHAPE+ARRAY=RESHAPE 6     B      MAPCUT_DBLE%TRANSPOSE+ARRAY=TRANSPOSE '   Ø     a   MAPCUT_DBLE%IMAP+ARRAY %   d  ¤   a   MAPCUT_DBLE%MM+ARRAY '        a   MAPCUT_DBLE%OMAP+ARRAY %     ¤   a   MAPCUT_DBLE%NN+ARRAY #   8         MAPCUT_CMPLX+ARRAY -   Ó  =      MAPCUT_CMPLX%DBLE+ARRAY=DBLE /     >      MAPCUT_CMPLX%AIMAG+ARRAY=AIMAG (   N     a   MAPCUT_CMPLX%IMAP+ARRAY &   Ú  ¤   a   MAPCUT_CMPLX%MM+ARRAY (   ~     a   MAPCUT_CMPLX%OMAP+ARRAY &   
  ¤   a   MAPCUT_CMPLX%NN+ARRAY ,   ®  ¬       SPIN_WEIGHT_1DARRAY+ANAFLAT 4   Z  <      SPIN_WEIGHT_1DARRAY%ABS+ANAFLAT=ABS 6     =      SPIN_WEIGHT_1DARRAY%DBLE+ANAFLAT=DBLE /   Ó     a   SPIN_WEIGHT_1DARRAY%SW+ANAFLAT /   _     a   SPIN_WEIGHT_1DARRAY%NN+ANAFLAT .   ó     a   SPIN_WEIGHT_1DARRAY%D+ANAFLAT 2     @   a   SPIN_WEIGHT_1DARRAY%TRANS+ANAFLAT .   Ç  @   a   SPIN_WEIGHT_1DARRAY%S+ANAFLAT ,     ¬       SPIN_WEIGHT_2DARRAY+ANAFLAT 4   ³  <      SPIN_WEIGHT_2DARRAY%ABS+ANAFLAT=ABS 6   ï  =      SPIN_WEIGHT_2DARRAY%DBLE+ANAFLAT=DBLE /   ,  ¤   a   SPIN_WEIGHT_2DARRAY%SW+ANAFLAT /   Ð     a   SPIN_WEIGHT_2DARRAY%NN+ANAFLAT .   d     a   SPIN_WEIGHT_2DARRAY%D+ANAFLAT 2   ø  @   a   SPIN_WEIGHT_2DARRAY%TRANS+ANAFLAT .   8  @   a   SPIN_WEIGHT_2DARRAY%S+ANAFLAT -   x        GAUSSIAN_ALM_1DARRAY+ANAFLAT =     @      GAUSSIAN_ALM_1DARRAY%PRESENT+ANAFLAT=PRESENT 5   K  <      GAUSSIAN_ALM_1DARRAY%ABS+ANAFLAT=ABS 5     <      GAUSSIAN_ALM_1DARRAY%INT+ANAFLAT=INT 7   Ã  =      GAUSSIAN_ALM_1DARRAY%DBLE+ANAFLAT=DBLE 9      >      GAUSSIAN_ALM_1DARRAY%DSQRT+ANAFLAT=DSQRT 7   >  =      GAUSSIAN_ALM_1DARRAY%SIZE+ANAFLAT=SIZE 5   {  <      GAUSSIAN_ALM_1DARRAY%MOD+ANAFLAT=MOD 9   ·  >      GAUSSIAN_ALM_1DARRAY%CMPLX+ANAFLAT=CMPLX 9   õ  >      GAUSSIAN_ALM_1DARRAY%CONJG+ANAFLAT=CONJG 0   3     a   GAUSSIAN_ALM_1DARRAY%NN+ANAFLAT /   Ç     a   GAUSSIAN_ALM_1DARRAY%D+ANAFLAT 0   [     a   GAUSSIAN_ALM_1DARRAY%IL+ANAFLAT 1   ï     a   GAUSSIAN_ALM_1DARRAY%ALM+ANAFLAT 0   {      a   GAUSSIAN_ALM_1DARRAY%CL+ANAFLAT 1   !  @   a   GAUSSIAN_ALM_1DARRAY%FIX+ANAFLAT -   G!        GAUSSIAN_ALM_2DARRAY+ANAFLAT =   Ú"  @      GAUSSIAN_ALM_2DARRAY%PRESENT+ANAFLAT=PRESENT 5   #  <      GAUSSIAN_ALM_2DARRAY%ABS+ANAFLAT=ABS 5   V#  <      GAUSSIAN_ALM_2DARRAY%INT+ANAFLAT=INT 7   #  =      GAUSSIAN_ALM_2DARRAY%DBLE+ANAFLAT=DBLE 9   Ï#  >      GAUSSIAN_ALM_2DARRAY%DSQRT+ANAFLAT=DSQRT 7   $  =      GAUSSIAN_ALM_2DARRAY%SIZE+ANAFLAT=SIZE 5   J$  <      GAUSSIAN_ALM_2DARRAY%MOD+ANAFLAT=MOD 9   $  >      GAUSSIAN_ALM_2DARRAY%CMPLX+ANAFLAT=CMPLX 9   Ä$  >      GAUSSIAN_ALM_2DARRAY%CONJG+ANAFLAT=CONJG 0   %     a   GAUSSIAN_ALM_2DARRAY%NN+ANAFLAT /   %     a   GAUSSIAN_ALM_2DARRAY%D+ANAFLAT 0   *&     a   GAUSSIAN_ALM_2DARRAY%IL+ANAFLAT 1   ¾&  ¤   a   GAUSSIAN_ALM_2DARRAY%ALM+ANAFLAT 0   b'     a   GAUSSIAN_ALM_2DARRAY%CL+ANAFLAT 1   î'  @   a   GAUSSIAN_ALM_2DARRAY%FIX+ANAFLAT    .(  «       DFT_1DARRAY     Ù(  <      DFT_1DARRAY%ABS     )  <      DFT_1DARRAY%SUM !   Q)  =      DFT_1DARRAY%DBLE     )     a   DFT_1DARRAY%MAP    *  ¤   a   DFT_1DARRAY%NN    ¾*  ¤   a   DFT_1DARRAY%D "   b+  @   a   DFT_1DARRAY%TRANS    ¢+  «       DFT_2DARRAY     M,  <      DFT_2DARRAY%ABS     ,  <      DFT_2DARRAY%SUM !   Å,  =      DFT_2DARRAY%DBLE     -  ¤   a   DFT_2DARRAY%MAP    ¦-     a   DFT_2DARRAY%NN    :.     a   DFT_2DARRAY%D "   Î.  @   a   DFT_2DARRAY%TRANS     /  Â       DFT_POL_1DARRAY &   Ð/  >      DFT_POL_1DARRAY%CONJG &   0  >      DFT_POL_1DARRAY%AIMAG %   L0  =      DFT_POL_1DARRAY%DBLE #   0  ¤   a   DFT_POL_1DARRAY%QU #   -1     a   DFT_POL_1DARRAY%NN "   Á1     a   DFT_POL_1DARRAY%D #   U2  ¤   a   DFT_POL_1DARRAY%EB &   ù2  @   a   DFT_POL_1DARRAY%TRANS     93  Â       DFT_POL_2DARRAY &   û3  >      DFT_POL_2DARRAY%CONJG &   94  >      DFT_POL_2DARRAY%AIMAG %   w4  =      DFT_POL_2DARRAY%DBLE #   ´4  ¼   a   DFT_POL_2DARRAY%QU #   p5     a   DFT_POL_2DARRAY%NN "   6     a   DFT_POL_2DARRAY%D #   6  ¼   a   DFT_POL_2DARRAY%EB &   T7  @   a   DFT_POL_2DARRAY%TRANS     7  Ä       DFT_ALL_1DARRAY &   X8  >      DFT_ALL_1DARRAY%CONJG &   8  >      DFT_ALL_1DARRAY%AIMAG %   Ô8  =      DFT_ALL_1DARRAY%DBLE $   9  ¤   a   DFT_ALL_1DARRAY%TQU #   µ9     a   DFT_ALL_1DARRAY%NN "   I:     a   DFT_ALL_1DARRAY%D $   Ý:  ¤   a   DFT_ALL_1DARRAY%TEB &   ;  @   a   DFT_ALL_1DARRAY%TRANS     Á;  Ä       DFT_ALL_2DARRAY &   <  >      DFT_ALL_2DARRAY%CONJG &   Ã<  >      DFT_ALL_2DARRAY%AIMAG %   =  =      DFT_ALL_2DARRAY%DBLE $   >=  ¼   a   DFT_ALL_2DARRAY%TQU #   ú=     a   DFT_ALL_2DARRAY%NN "   >     a   DFT_ALL_2DARRAY%D $   "?  ¼   a   DFT_ALL_2DARRAY%TEB &   Þ?  @   a   DFT_ALL_2DARRAY%TRANS    @         DERIVEMAP_ALL #   ´@  =      DERIVEMAP_ALL%SIZE #   ñ@  =      DERIVEMAP_ALL%DBLE     .A     a   DERIVEMAP_ALL%D "   ÂA  ¤   a   DERIVEMAP_ALL%MAP '   fB  Ô   a   DERIVEMAP_ALL%DERIVMAP !   :C  ª       DERIVEMAP_NTH_1D '   äC  >      DERIVEMAP_NTH_1D%CMPLX &   "D  =      DERIVEMAP_NTH_1D%DBLE $   _D     a   DERIVEMAP_NTH_1D%NN #   óD     a   DERIVEMAP_NTH_1D%D %   E     a   DERIVEMAP_NTH_1D%MAP &   F  ¤   a   DERIVEMAP_NTH_1D%DMAP %   ·F  @   a   DERIVEMAP_NTH_1D%NTH !   ÷F  ¥       DERIVEMAP_NTH_2D &   G  =      DERIVEMAP_NTH_2D%SIZE &   ÙG  =      DERIVEMAP_NTH_2D%DBLE #   H     a   DERIVEMAP_NTH_2D%D %   ªH  ¤   a   DERIVEMAP_NTH_2D%MAP *   NI  ¼   a   DERIVEMAP_NTH_2D%DERIVMAP %   
J  @   a   DERIVEMAP_NTH_2D%NTH %   JJ         DERIVEMAP_ALM_NTH_1D *   ÜJ  =      DERIVEMAP_ALM_NTH_1D%DBLE (   K     a   DERIVEMAP_ALM_NTH_1D%NN '   ­K     a   DERIVEMAP_ALM_NTH_1D%D )   AL     a   DERIVEMAP_ALM_NTH_1D%ALM *   ÍL  ¤   a   DERIVEMAP_ALM_NTH_1D%DMAP )   qM  @   a   DERIVEMAP_ALM_NTH_1D%NTH    ±M         QU2EB    3N  =      QU2EB%DBLE    pN     a   QU2EB%NN    O     a   QU2EB%D    O  ¤   a   QU2EB%QU    <P  ¤   a   QU2EB%EB    àP  @   a   QU2EB%TRANS     Q         PUREEB    »Q  @      PUREEB%PRESENT    ûQ  =      PUREEB%DBLE    8R  ¤   a   PUREEB%QU    ÜR  ¤   a   PUREEB%NN    S  ¤   a   PUREEB%D    $T  ¤   a   PUREEB%EB    ÈT     a   PUREEB%W    TU  ¤   a   PUREEB%WD    øU         ARRAY_DERIV !   V  =      ARRAY_DERIV%DBLE    ÔV     a   ARRAY_DERIV%NN    hW     a   ARRAY_DERIV%D    üW     a   ARRAY_DERIV%W    X     a   ARRAY_DERIV%WX    Y     a   ARRAY_DERIV%WY      Y     a   ARRAY_DERIV%WXX     ,Z     a   ARRAY_DERIV%WXY     ¸Z     a   ARRAY_DERIV%WYY    D[         GAUSSIAN_MAP %   ß[  @      GAUSSIAN_MAP%PRESENT     \  ¤   a   GAUSSIAN_MAP%NN     Ã\  ¤   a   GAUSSIAN_MAP%MM    g]  ¤   a   GAUSSIAN_MAP%D     ^  ¤   a   GAUSSIAN_MAP%EL     ¯^     a   GAUSSIAN_MAP%CL !   ;_     a   GAUSSIAN_MAP%MAP !   Ç_  @   a   GAUSSIAN_MAP%FIX !   `  ­       GAUSSIAN_MAP_POL )   ´`  @      GAUSSIAN_MAP_POL%PRESENT $   ô`  ¤   a   GAUSSIAN_MAP_POL%NN $   a  ¤   a   GAUSSIAN_MAP_POL%MM #   <b  ¤   a   GAUSSIAN_MAP_POL%D $   àb  ¤   a   GAUSSIAN_MAP_POL%EL $   c     a   GAUSSIAN_MAP_POL%EE $   d     a   GAUSSIAN_MAP_POL%BB $   d  ¤   a   GAUSSIAN_MAP_POL%QU #   @e     a   GAUSSIAN_MAP_POL%P %   Ìe  @   a   GAUSSIAN_MAP_POL%FIX 