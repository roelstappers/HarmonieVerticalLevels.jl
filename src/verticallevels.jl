using OffsetArrays

# JPN: Total number of layers.
JPN=55
JPNP=JPN+1
# JPNPBL: Total number of layers in the PBL.
 JPNPBL=8
# JPNSTR: Total number of layers in the stratosphere.
JPNSTR=27
# JPNPRES: Total number of pure pressure layers (minimum 1 included).
 JPNPRES=12
!#JPNSIGM: Total number of pure sigma layers (minimum 1 included).
JPNSIGM=1

ZM= OffsetArray{Float64}(undef,0:JPN)
ZH= OffsetArray{Float64}(undef,0:JPN)
ZM .= 0.0
ZH .= 0.0
ZAF = zeros(JPNP)
ZBF = zeros(JPNP)
ZRST  = zeros(JPNP)
ZSIG = zeros(JPNP)

IN=JPN
INP=JPNP
INPBL=JPNPBL
INSTR=JPNSTR
INPRES=JPNPRES
INSIGM=JPNSIGM 

ZIN=IN
ZEPSD = 0.0001


#     * ZPRE_PBL: pressure of the top of the PBL in pascals:
ZPRE_PBL=90000.
#     * ZPRE_STR: pressure of the tropopause in pascals:
ZPRE_STR=25000.
# ZP1: pressure of the full layer l=1, in pascals:
ZP1=9.9
#  ZDPRE_FLAYL: [Delta pressure] of the full layer l=L, in pascals:
ZDPRE_FLAYL=410.
#   * ZVP00: standard surface pressure, in pascals:
ZVP00=101325.
# ZPREF: standard pressure at 200m (mean orography on the Earth),
#       in pascals:
ZPREF=98945.37974

#!     * LLAPRXPK:
#!       Full layers are assumed to be computed as for the options
#!       LVERTFE=F, NDLNPR=0 of ARPEGE/ALADIN.
#!       LLAPRXPK=T => pressure(l)=0.5(pressure(lbar-1)+pressure(lbar))
#!        ("l" stands for full levels, "lbar" for half levels).
#!       LLAPRXPK=F => a more tricky way to compute pressure(l).
#!       When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=F
#!        in the model, it is recommended to use LLAPRXPK=F.
#!       When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=T
#!        of for LVERTFE=T, it is recommended to use LLAPRXPK=T.
LLAPRXPK= true 
if LLAPRXPK
    ZP1H=2*ZP1
else
    ZP1H=exp(1.0)*ZP1
end 
    

#  ZALP1,ZALP3,ZALPH: coefficients defining the function allowing
#  to compute the A and B.
#  For ZALP1 and ZALP3 it is recommended to take values between 1 and 5.
#  For ZALPH it is recommended to take values between -3 and -1
# (ZALPH must bever be > -1).

ZALP1=2.8
ZALP3=1.7
ZALPH=-1.6 



ZX1=1.0/ZIN
ZX2=INSTR/ZIN
ZX3=(IN-INPBL)/ZIN
ZX4=(IN-1.0)/ZIN

# !     * ZXP3 to ZXP4:
ZXP3=1.0-ZX3
ZXP4=1.0-ZX4

#    * ZY1 to ZY4
ZY1=ZP1H/ZVP00
ZY2=ZPRE_STR/ZVP00 
ZY3=ZPRE_PBL/ZVP00 
ZY4=(ZVP00-ZDPRE_FLAYL)/ZVP00

# !     * ZYP3 to ZYP4:
ZYP3=1-ZY3
ZYP4=1-ZY4

#     =======================================================================
#     Computation of mapping function ZM(JLEV)=m(x(JLEV))

#     * ZM between the half-level 1 and the tropopause level:

ZDFAC1=(ZX1*ZY2-ZX2*ZY1)*(1/ZX1)*(ZX2-ZX1)^(-ZALP1)


for JLEV=1:INSTR
  ZX=JLEV/ZIN
  ZM[JLEV]=ZX*ZY1/ZX1+ZDFAC1*(ZX-ZX1)^ZALP1
end 

#     * ZM derivative at tropopause level:

ZX2M  = INSTR/ZIN  - ZEPSD
ZX2MM = INSTR/ZIN - 2 * ZEPSD
ZM2M  = ZX2M*ZY1/ZX1+ZDFAC1*(ZX2M-ZX1)^ZALP1
ZM2MM = ZX2MM*ZY1/ZX1+ZDFAC1*(ZX2MM-ZX1)^ZALP1

ZDERI1=(ZM[INSTR]-ZM2M)/ZEPSD
ZDERI11=(ZM[INSTR]-2*ZM2M+ZM2MM)/(ZEPSD*ZEPSD)

#     * ZM between the PBL top and the half-level IN.

ZAA1=(((0.8*ZY3-ZY2)/(ZX3-ZX2))-(ZY1/ZX1)) *(ZX1*(ZX2-ZX1)/(ZX1*ZY2-ZX2*ZY1))
ZAA3=(((1.4*ZY2-ZY3)/(ZX2-ZX3))-(ZYP4/ZXP4)) *(ZXP4*(ZXP3-ZXP4)/(ZXP4*ZYP3-ZXP3*ZYP4))

println(" Proposed alpha1: $ZAA1")
println(" Proposed alpha3: $ZAA3")

ZDFAC3=(ZXP4*ZYP3-ZXP3*ZYP4)*(1/ZXP4)*(ZX4-ZX3)^(-ZALP3)

#m(x) = 1.0-(1.0-x)*(ZYP4/ZXP4)-ZDFAC3*(ZX4-x)^ZALP3
#JLEV = (IN-INPBL):(IN-1)
# ZM = m.(JLEV./ZIN)

for  JLEV = (IN-INPBL):(IN-1)
    ZX=JLEV/ZIN
    ZM[JLEV]=1-(1-ZX)*(ZYP4/ZXP4)-ZDFAC3*(ZX4-ZX)^ZALP3

end 

#     * ZM derivative at PBL top:

      ZX3P= (IN-INPBL)/ZIN + ZEPSD
      ZX3PP= (IN-INPBL)/ZIN + 2* ZEPSD
      ZM3P = 1-(1-ZX3P)*(ZYP4/ZXP4)-ZDFAC3*(ZX4-ZX3P)^ZALP3
      ZM3PP = 1-(1-ZX3PP)*(ZYP4/ZXP4)-ZDFAC3*(ZX4-ZX3PP)^ZALP3

      ZDERI2=(ZM3P-ZM[IN-INPBL])/ZEPSD
      ZDERI22=(ZM3PP-2*ZM3P +ZM[IN-INPBL])/(ZEPSD*ZEPSD)

#     * ZM between the stratosphere level and the PBL top.

      ZDX=ZX3-ZX2
      ZDY=ZY3-ZY2
      ZS=ZDY/ZDX

      for JLEV=INSTR+1:(IN-INPBL)-1
        ZX=JLEV/ZIN
        ZM[JLEV]=ZY2+(ZX-ZX2)*ZDERI1   +(ZDX*(ZS-ZDERI1)+(ZX-ZX3)*(ZDERI1+ZDERI2-2*ZS)) *(ZX-ZX2)*(ZX-ZX2)/(ZDX*ZDX)
      end 

#     * ZM between the top and the half-level 1:

      for JLEV=0:1
        ZX=JLEV/ZIN
        ZM[JLEV]=ZX*ZY1/ZX1
      end 

#     * ZM between the half-level IN-1 and the surface:

      for JLEV=(IN-1):IN
        ZX=JLEV/ZIN
        ZM[JLEV]=1-(1-ZX)*((1-ZY4)/(1-ZX4))
      end 

#     * Bound ZM between 0 and 1:

      for  JLEV=0:IN
        ZM[JLEV]=maximum([0.,minimum([1.0,ZM[JLEV]])])
      end 

      println(" ZM: ")
      for JLEV=0:IN
        println(ZM[JLEV])
      end 


#     =======================================================================
#     Computation of *mapped hybridicity* function ZH(JLEV)=h(m(x(JLEV)))

      ZETAP=ZM[INPRES]
      ZETAS=ZM[IN-INSIGM]

      ZAA=ZALPH*ZETAS*ZETAS/(ZETAS-ZETAP)
      ZBB=1+ZAA/ZETAS

      println(" ZETAP: $ZETAP")
      println(" ZETAS: $ZETAS")
      println(" ZAA: $ZAA")
      println(" ZBB: $ZBB")

      for  JLEV=0:INPRES
        ZH(JLEV)=0.0
      end

      for JLEV=INPRES+1:(IN-INSIGM)-1
        ZX=ZM[JLEV]
        ZH[JLEV]=ZAA/(ZBB-((ZX-ZETAP)/(ZETAS-ZETAP))^ZALPH)
      end 

      for JLEV=(IN-INSIGM):IN
        ZX=ZM[JLEV]
        ZH[JLEV]=ZX
      end 

      for JLEV=0:IN
        ZH[JLEV]=maximum([0.,minimum([1.0,ZH[JLEV]])])
      end 

      println(" ZH: ")
      for JLEV=0:IN
        println(ZH[JLEV])
      end 

#     =======================================================================

#     * A and B functions on half layers (put in ZAF and ZBF):

      for  JLEV=0:IN
        ZX=JLEV/ZIN
        ZAF[JLEV+1]=ZVP00*(ZM[JLEV]-ZH[JLEV])
        ZBF[JLEV+1]=ZH[JLEV]
      end 

#     =======================================================================

#     * Printings:

      for JK=1:IN+1
        ZRST[JK]=ZBF[JK]*ZPREF/(ZBF[JK]*ZPREF+ZAF[JK])
        ZSIG[JK]=ZAF[JK]/ZPREF+ZBF[JK]
      end 

      
      println("Number of levels= $IN") 
      println("Number of hybrid levels= $(IN-INPRES-INSIGM)")
      println("Number of pure sigma levels= $INSIGM")
      println("Number of pure pressure levels= $INPRES")
      println("Number of levels in the stratosphere= $INSTR")
      println("Number of levels in the PBL=$INPBL")
      println("Reference pressure at 200m: ZPREF= $ZPREF")
      println("Reference pressure at 0m: ZVP00=$ZVP00")
      println("Reference pressure of the top of the PBL: ZPRE_PBL= $ZPRE_PBL")
      println("Reference pressure of the tropopause: ZPRE_STR= $ZPRE_STR")
      println("Pressure of the layer nr 1: ZP1=$ZP1")
      println("Depth pressure of the layer nr $IN")
      println("ZDPRE_FLAYL=$ZDPRE_FLAYL")
      println("LLAPRXPK=$LLAPRXPK")
      println("ZALP1=$ZALP1")
      println("ZALP3=$ZALP3")
      println("ZALPH=$ZALPH")
      println("ILa      A            B          Sigma      1 - Sigma   Rap Sig-Hyb Prehyd(lbar)    Prehyd(l)")

      for JK=1:1
#       println(""  JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),
#     &   ZAF(JK)+ZBF(JK)*ZVP00
#      ENDDO
#      IF (.NOT.LLAPRXPK) THEN
#        DO JK=2,2
#          WRITE(IFILE1,'(I3,F13.6,F13.10,3F13.7,F16.6,F16.6)') 
#     &     JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),
#     &     ZAF(JK)+ZBF(JK)*ZVP00,
#     &     (ZAF(JK)+ZBF(JK)*ZVP00)/EXP(1.)
#        ENDDO
#        DO JK=3,IN+1
#          ZPREHYDHM=ZAF(JK-1)+ZBF(JK-1)*ZVP00
#          ZPREHYDH =ZAF(JK)+ZBF(JK)*ZVP00
#          ZPREHYDF=EXP(
#     &    (ZPREHYDH*LOG(ZPREHYDH)-ZPREHYDHM*LOG(ZPREHYDHM))/
#     &    (ZPREHYDH-ZPREHYDHM) - 1
#     &    )
#          WRITE(IFILE1,'(I3,F13.6,F13.10,3F13.7,F16.6,F16.6)') 
#     &     JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),
#     &     ZAF(JK)+ZBF(JK)*ZVP00,ZPREHYDF
#        ENDDO
#      ELSE
#        DO JK=2,2
#          WRITE(IFILE1,'(I3,F13.6,F13.10,3F13.7,F16.6,F16.6)') 
#     &     JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),
#     &     ZAF(JK)+ZBF(JK)*ZVP00,
#     &     (ZAF(JK)+ZBF(JK)*ZVP00)/2.
#        ENDDO
#        DO JK=3,IN+1
#          ZPREHYDHM=ZAF(JK-1)+ZBF(JK-1)*ZVP00
#          ZPREHYDH =ZAF(JK)+ZBF(JK)*ZVP00
#          ZPREHYDF=0.5*(ZPREHYDHM+ZPREHYDH)
#          WRITE(IFILE1,'(I3,F13.6,F13.10,3F13.7,F16.6,F16.6)') 
#     &     JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),
#     &     ZAF(JK)+ZBF(JK)*ZVP00,ZPREHYDF
#        ENDDO
#      ENDIF
#
#      WRITE(IFILE2,*)
#      DO JK=1,IN+1
#        WRITE(IFILE2,'(A12,F9.5,A28,I2)')
#     &   '\\put( 20.00,',100.*(1.-ZSIG(JK)),
#     &   '){\\line(1,0){50.00}} % jlev=',JK-1
#      ENDDO
#
#      WRITE(IFILE3,'(A)') ' Namelist obtained with:'
#
#      WRITE(IFILE3,*)
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of levels=',IN
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of hybrid levels=',IN-INPRES-INSIGM
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of pure sigma levels=',INSIGM
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of pure pressure levels=',INPRES
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of levels in the stratosphere=',INSTR
#      WRITE(IFILE3,'(1X,A,1X,I4)')
#     & ' * Number of levels in the PBL=',INPBL
#
#      WRITE(IFILE3,*)
#      WRITE(IFILE3,'(1X,A,1X,F15.7)')
#     & ' * Reference pressure at 200m: ZPREF=',ZPREF
#      WRITE(IFILE3,'(1X,A,1X,F15.7)')
#     & ' * Reference pressure at 0m: ZVP00=',ZVP00
#      WRITE(IFILE3,'(1X,A,1X,F15.7)')
#     & ' * Reference pressure of the top of the PBL: ZPRE_PBL=',ZPRE_PBL
#      WRITE(IFILE3,'(1X,A,1X,F15.7)')
#     & ' * Reference pressure of the tropopause: ZPRE_STR=',ZPRE_STR
#      WRITE(IFILE3,'(1X,A,1X,F15.7)')
#     & ' * Pressure of the layer nr 1: ZP1=',ZP1
#      WRITE(IFILE3,'(1X,A,I3,A,1X,F15.7)')
#     & ' * Depth pressure of the layer nr ',IN,
#     & ': ZDPRE_FLAYL=',ZDPRE_FLAYL
#
#      WRITE(IFILE3,*)
#      WRITE(IFILE3,'(1X,A,1X,L2)') ' * LLAPRXPK=',LLAPRXPK
#      WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * ZALP1=',ZALP1
#      WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * ZALP3=',ZALP3
#      WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * ZALPH=',ZALPH
#
#      WRITE(IFILE3,*)
#      WRITE(IFILE3,*) ' &NAMVV1'
#      WRITE(IFILE3,*) '    DVALH(0)=0.,'
#      DO JK=2,IN+1
#        WRITE(IFILE3,'(5X,F12.6,A)') ZAF(JK),','
#      ENDDO
#      WRITE(IFILE3,*) '    DVBH(0)=0.,'
#      DO JK=2,IN+1
#        WRITE(IFILE3,'(5X,F12.10,A)') ZBF(JK),','
#      ENDDO
#