	  �%  V   k820309              12.0        ��AY                                                                                                           
       mylapack95.f90 MYLAPACK95                                                    
       #         @                                      	               #A    #IPIV    #RCOND    #NORM    #INFO              
                                                  	 �             &                   &                                                                                                         �             &                                                                                         	                 
                                                                                                                #         @                                     	               #A 	   #IPIV 
   #RCOND    #NORM    #INFO              
                               	                   
 �             &                   &                                                                                     
                    �             &                                                                                         
                 
                                                                                                                #         @                                      	               #A    #IPIV    #RCOND    #NORM    #INFO              
                                                   �             &                   &                                                                                                         �             &                                                                                         	                 
                                                                                                                #         @                                      	               #A    #IPIV    #RCOND    #NORM    #INFO              
                                                   �             &                   &                                                                                                         �             &                                                                                         
                 
                                                                                                                #         @                                      	               #A    #IPIV    #INFO              
                                                  	 �             &                   &                                                     
                                                     �            &                                                                                                 #         @                                     	               #A    #IPIV     #INFO !             
                                                  
 �             &                   &                                                     
                                                      �            &                                                                                     !            #         @                                  "    	               #A #   #IPIV $   #INFO %             
                               #                    �             &                   &                                                     
                                 $                    �            &                                                                                     %            #         @                                  &    	               #A '   #IPIV (   #INFO )             
                               '                    �             &                   &                                                     
                                 (                    �            &                                                                                     )            #         @                                  *    	               #A +   #B ,   #W -   #ITYPE .   #JOBZ /   #UPLO 0   #INFO 1             
                               +                   	 Q             &                   &                                                     
                               ,                   	 R             &                   &                                                                                    -                   	 S             &                                                     
                                .                     
                               /                                     
                               0                                                                     1            #         @                                 2    	               #A 3   #B 4   #W 5   #ITYPE 6   #JOBZ 7   #UPLO 8   #INFO 9             
                               3                   
 T             &                   &                                                     
                               4                   
 U             &                   &                                                                                    5                   
 V             &                                                     
                                6                     
                               7                                     
                               8                                                                     9            #         @                                  :                  #INV_LAPACK%PRESENT ;   #INV_LAPACK%SIZE <   #K =   #M >   #ACCURACY ?                                              ;     PRESENT                                            <     SIZE           
D@                               =                   
               &                   &                                                     F @                               >                   
               &                   &                                                     
 @                               ?     
      #         @                                 @                  #CHECK_INV_LAPACK%ABS A   #CHECK_INV_LAPACK%MATMUL B   #CHECK_INV_LAPACK%SIZE C   #M D   #IM E   #ACCURACY F                                              A     ABS                                            B     MATMUL                                            C     SIZE           
 @                               D                   
              &                   &                                                     
                                  E                   
              &                   &                                                     
                                  F     
      #         @                                  G                  #SYGV_LAPACK%PRESENT H   #SYGV_LAPACK%SIZE I   #K J   #U K   #D L   #M M                                              H     PRESENT                                            I     SIZE           
                                 J                   
 
              &                   &                                                     D @                               K                   
               &                   &                                                     D@                               L                   
               &                                                      @                               M                   
               &                   &                                           #         @                                  N                  #CHECK_SYGV_LAPACK%TRANSPOSE O   #CHECK_SYGV_LAPACK%ABS P   #CHECK_SYGV_LAPACK%MATMUL Q   #CHECK_SYGV_LAPACK%SIZE R   #D S   #U T   #M U                                              O     TRANSPOSE                                            P     ABS                                            Q     MATMUL                                            R     SIZE           
                                  S                   
              &                                                     
  @                               T                   
              &                   &                                                     
 @                               U                   
              &                   &                                              �   "      fn#fn    �   @   J   F95_LAPACK &     x       SGETRF_F95+F95_LAPACK (   z  �   a   SGETRF_F95%A+F95_LAPACK +     �   a   SGETRF_F95%IPIV+F95_LAPACK ,   �  @   a   SGETRF_F95%RCOND+F95_LAPACK +   �  P   a   SGETRF_F95%NORM+F95_LAPACK +   :  @   a   SGETRF_F95%INFO+F95_LAPACK &   z  x       DGETRF_F95+F95_LAPACK (   �  �   a   DGETRF_F95%A+F95_LAPACK +   �  �   a   DGETRF_F95%IPIV+F95_LAPACK ,   "  @   a   DGETRF_F95%RCOND+F95_LAPACK +   b  P   a   DGETRF_F95%NORM+F95_LAPACK +   �  @   a   DGETRF_F95%INFO+F95_LAPACK &   �  x       CGETRF_F95+F95_LAPACK (   j  �   a   CGETRF_F95%A+F95_LAPACK +     �   a   CGETRF_F95%IPIV+F95_LAPACK ,   �  @   a   CGETRF_F95%RCOND+F95_LAPACK +   �  P   a   CGETRF_F95%NORM+F95_LAPACK +   *  @   a   CGETRF_F95%INFO+F95_LAPACK &   j  x       ZGETRF_F95+F95_LAPACK (   �  �   a   ZGETRF_F95%A+F95_LAPACK +   �	  �   a   ZGETRF_F95%IPIV+F95_LAPACK ,   
  @   a   ZGETRF_F95%RCOND+F95_LAPACK +   R
  P   a   ZGETRF_F95%NORM+F95_LAPACK +   �
  @   a   ZGETRF_F95%INFO+F95_LAPACK &   �
  c       SGETRI_F95+F95_LAPACK (   E  �   a   SGETRI_F95%A+F95_LAPACK +   �  �   a   SGETRI_F95%IPIV+F95_LAPACK +   u  @   a   SGETRI_F95%INFO+F95_LAPACK &   �  c       DGETRI_F95+F95_LAPACK (     �   a   DGETRI_F95%A+F95_LAPACK +   �  �   a   DGETRI_F95%IPIV+F95_LAPACK +   H  @   a   DGETRI_F95%INFO+F95_LAPACK &   �  c       CGETRI_F95+F95_LAPACK (   �  �   a   CGETRI_F95%A+F95_LAPACK +   �  �   a   CGETRI_F95%IPIV+F95_LAPACK +     @   a   CGETRI_F95%INFO+F95_LAPACK &   [  c       ZGETRI_F95+F95_LAPACK (   �  �   a   ZGETRI_F95%A+F95_LAPACK +   b  �   a   ZGETRI_F95%IPIV+F95_LAPACK +   �  @   a   ZGETRI_F95%INFO+F95_LAPACK %   .  �       SSYGV_F95+F95_LAPACK '   �  �   a   SSYGV_F95%A+F95_LAPACK '   X  �   a   SSYGV_F95%B+F95_LAPACK '   �  �   a   SSYGV_F95%W+F95_LAPACK +   �  @   a   SSYGV_F95%ITYPE+F95_LAPACK *   �  P   a   SSYGV_F95%JOBZ+F95_LAPACK *     P   a   SSYGV_F95%UPLO+F95_LAPACK *   h  @   a   SSYGV_F95%INFO+F95_LAPACK %   �  �       DSYGV_F95+F95_LAPACK '   .  �   a   DSYGV_F95%A+F95_LAPACK '   �  �   a   DSYGV_F95%B+F95_LAPACK '   v  �   a   DSYGV_F95%W+F95_LAPACK +     @   a   DSYGV_F95%ITYPE+F95_LAPACK *   B  P   a   DSYGV_F95%JOBZ+F95_LAPACK *   �  P   a   DSYGV_F95%UPLO+F95_LAPACK *   �  @   a   DSYGV_F95%INFO+F95_LAPACK    "  �       INV_LAPACK #   �  @      INV_LAPACK%PRESENT     �  =      INV_LAPACK%SIZE    0  �   a   INV_LAPACK%K    �  �   a   INV_LAPACK%M $   x  @   a   INV_LAPACK%ACCURACY !   �  �       CHECK_INV_LAPACK %   o  <      CHECK_INV_LAPACK%ABS (   �  ?      CHECK_INV_LAPACK%MATMUL &   �  =      CHECK_INV_LAPACK%SIZE #   '  �   a   CHECK_INV_LAPACK%M $   �  �   a   CHECK_INV_LAPACK%IM *   o  @   a   CHECK_INV_LAPACK%ACCURACY    �  �       SYGV_LAPACK $   B  @      SYGV_LAPACK%PRESENT !   �  =      SYGV_LAPACK%SIZE    �  �   a   SYGV_LAPACK%K    c   �   a   SYGV_LAPACK%U    !  �   a   SYGV_LAPACK%D    �!  �   a   SYGV_LAPACK%M "   7"  �       CHECK_SYGV_LAPACK ,   
#  B      CHECK_SYGV_LAPACK%TRANSPOSE &   L#  <      CHECK_SYGV_LAPACK%ABS )   �#  ?      CHECK_SYGV_LAPACK%MATMUL '   �#  =      CHECK_SYGV_LAPACK%SIZE $   $  �   a   CHECK_SYGV_LAPACK%D $   �$  �   a   CHECK_SYGV_LAPACK%U $   4%  �   a   CHECK_SYGV_LAPACK%M 