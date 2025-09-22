!========================================================
program Double_transfer_v1
implicit none
!All parameters passed through from input file
!Modified Lax - adjusted weighting on central node.  Equal weighting is too dissipative
!Finds maximum porosity and reports time dependency of T, Cb, Cl, Ib, Il at this location
!Work on convection parked in this version.  Includes instead:
!Nd 143 versus Nd 144 [done]
!Melt viscosity depends on composition [done]
!Under versus overaccretion [done]
!Modified matrix viscosity [implemented here]
!Modified random number approach compared to version 16
!Modifed emplacement approach compared to version 16
!Density and/or porosity controlled emplacement with random element
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Counters
integer i,j,k
integer n
parameter (n=1000000)
integer velflag
real(8) runtime

!Arrays
real(8) ha(-n:n),hb(-n:n)	    !Enthalpy
real(8) ta(-n:n),tb(-n:n)	    !Temperature
real(8) pa(-n:n),pb(-n:n)	    !Porosity
real(8) ca(-n:n),cb(-n:n)	    !Bulk Composition
real(8) cc(-n:n)                !Bulk composition for checking
real(8) ws(-n:n),wl(-n:n)	    !Velocity
real(8) csa(-n:n),csb(-n:n)	    !Solid Composition
real(8) cla(-n:n),clb(-n:n)	    !Liquid Composition
real(8) ica(-n:n),icb(-n:n)	    !Isotope Bulk Composition 1
real(8) icc(-n:n)       	    !Isotope Bulk Composition 1 for checking
real(8) icsa(-n:n),icsb(-n:n)	!Isotope Solid Composition 1
real(8) icla(-n:n),iclb(-n:n)	!Isotope Liquid Composition 1
real(8) ima(-n:n),imb(-n:n)	    !Isotope Bulk Composition 2
real(8) imsa(-n:n),imsb(-n:n)	!Isotope Solid Composition 2
real(8) imla(-n:n),imlb(-n:n)	!Isotope Liquid Composition 2
real(8) imc(-n:n)               !isotope check 2
real(8) melt(-n:n)			    !Melting Rate
real(8) clostim(-n:n)           !Closure time
real(8) density(-n:n)           !Density
real(8) denscont(-n:n)          !Density contrast        
!integer layerbase(1:100)
!integer layertop(1:100)
real(8) htotal(1:100)
real(8) have(1:100)
real(8) ctotal(1:100)
real(8) cave(1:100)
real(8) cltotal(1:100)
real(8) clave(1:100)
real(8) cstotal(1:100)
real(8) csave(1:100)
real(8) crysum(1:100)
integer advecttype(-n:n)
integer icum(1:100)
real(8) dcum(1:100)
real(8) Melt_seg(-n:n)
real(8) pmag(-n:n) 

!Functions
real(8) enth_dim_to_nondim,temp_dim_to_nondim
real(8) enth_nondim_to_dim,temp_nondim_to_dim
real(8) SolidComp,LiquidComp
real(8) LFraction
real(8) IsoSolid,IsoLiquid
real(8) Liquidus
real(8) ndsio

!File Input/Output
integer outcount,freqcount,timeoutfreq,timeoutcount
character*4 filecount
character*5 filetime
character*50 dummy
character*14 sillstyle
real(8) OutRate

!Real Parameters
real(8) a,b,alpha,beta,drho,rho,g,SLT,max_melt_viscosity,min_melt_viscosity,shear_viscosity,cee
real(8) kt,Lf,Cp,Tl,Hs,Hl,Ae,PartitionA,PartitionB
real(8) TSill,TSill_dummy,Sill_Width,Sill_Comp,Base,Top,Geotherm,Model_Time
real(8) Emplacement_Depth1,Emplacement_Width1,Emplacement_Rate1,dzn
real(8) Emplacement_Depth2,Emplacement_Width2,Emplacement_Rate2
real(8) Emplacement_Depth3,Emplacement_Width3,Emplacement_Rate3
real(8) A2,B2,C2,laxc,laxi,CloseTemp,silltime,dtsill,Output_Interval,dum
real(8) aSiO2,bSiO2,cSiO2,dSiO2,depthoutput,DensLiqB,DensLiqA,DensSolB,DensSolA,basdens

!ND Parameters
real(8) perm,delta,tau,omega,keff,stefan,Year,keff_dummy
real(8) Time,ND_Geotherm,ddtsill,dsilltime,temptop,weight,weighttemp 

integer z1,z2,layernumb,loop,depthi,testflag,empexp,multdepth
integer d11,d12,d21,d22,d13,d23,emplacementflagp,emplacementflagd,emplacementpi,emplacementdi
integer errflag

real(8) dtdiff,dtadv,dz,dt1,dt2,tmult,dtl,csum,isum,maxvel,dtartdiff,maxdiff,casum,ccsum,cerror,cumerror,ierror, imerror

!Initial Conditions
real(8) InitialH,InitialT,InitialP,CbU,CbL,InitialI1,InitialI2
!real(8) InitialCIS,InitialCIL

!Sill Paramters
integer SillNodes,SillN1,SillN2,SillN3,Accrete,SillCount1,SillCount2,randdepthn1,randdepthn2,SillRandom1,SillCount3,randdepthn3
integer initsill,SillRandom2,SillRandom3
real(8) SillEnth,SillTemp,SillPhi,SillComp,SillIso1,SillIso2,SillRate1,SillRate2,Randdepth1,Randdepth2,SillRate3,Randdepth3
real(8) Emplacement_Pause1, StartTime, Emplacement_Pause2

!Crust Parameters
real(8) Crust_Depth
integer CrustND

integer seed,d2rand,maxpi,maxti,convflag,sample,emplacementi,maxvi,maxdi
real(8) maxp,maxt,writetime,porodepth,diffmax,diffmin,diffexp,tmultdiff,tmultad1,tmultad2
integer random,densflag

!Transfer Parameters
integer j1, j2, j3, j4, volFlag, DMFn, Voln, Tdepth, Deni, outcountSM, IorErF,MMorMFlag
real(8) cpa, pressureE, volCh, volT, radg, dco
real(8) AvpaStart, OvCa, AvCa, OvH, AvH, AvCs, AvCl, Avdensill, AvPa, AvT
real(8) mpc, mc, t0, hmelt, D, PI,denR, hmushm
real(8) AREAU, AREAL, DMF
character*5 filecountSM, filecountsm1, FILECOUNTSM1L

!compsolve
integer baseflagU, baseflagL, topflagU, topflagL
real(8) cerrorU, cerrorL
real(8) outcomp, outcc

!Porosity
integer maxlpi, maxupi
real(8) maxlp, maxup

! Added 
integer TopFlagP, BaseFlagP, MaxPor, SN, jemh, j3new
integer xaboven
real EE, EmHMF, ovdenR, xabove

!! added 26/08/20
integer transon, transcount, transn, transnodes, e1, e2, dene, ovtn, outcountsm1
real(8) transwidth,Tran_rate,transrate, runtimeT, eavca,eavh,eavpa,eavcs,eavcl,eavdensill,eavt, occ, oh, oca, avcd

! added 28/08/20
integer e3, e4, e3new
real(8) epressuree, eovdenr, edenr, et0

! added 25/09
integer j4new, writeFL, TransF, e4new
real(8) ovdensmagL, ovclamagl, densmagl, clamagl, thermdiffl, viscmagl, dtwrl, JPBeta, DPL, elsc
real(8) ovdensmagu, ovclamagu, densmagu, clamagu, thermdiffu, viscmagu, dtwru, dpu
real(8), parameter::PI_8=4*atan(1.0_8)

!! added 02/10/20
real(8) t2,taue,tempdx,propl,hdot,hlm1,TotPressure,t1,lambda,grhl,HRTI0,HRTI,pressureRTI,hlm,hmb
real(8) crustmu
integer tempdxn,j3b,j4b,j1b,j2b,ff,transk

!! added 30/10/20
integer AmountS, AmountI, AmountE, AELower, AILower
real(8)  Uef

!!added 04/03/21
real(8) ndcheck, sio2check, cacheck, densscheck, denslcheck, visccheck

!added 9/03/21
real(8) lambdap, ovdensb, ovporb, avporb, MushShear, alphamei, Magvisc, avdensb

!added 10/03
real(8) ocab, ohb, avcab, avhb, avcsb, avclb, viscmagb, checkMEI
real(8) hrtiP

!added 3/6 to do with phase diagram
real(8) Ts(-n:n), A1(-n:n), B1(-n:n), C1(-n:n)
real(8) TsL, TsU, A1U, B1U, C1U, A1L, B1L, C1L
integer duln, TrnsInProg, L1,L2,L3,L4
real(8) dul,dcoL

!transfer parameters v2
real(8) pressureEL, totpressureL, tauEL, LDPL, t1L, t2L, HRTIPL, HRTIL, hlm1L, lambdaL
real(8) LambdaPL, OVDENRL, DENRL, HMBL, HLML, GRHLL, HDOTL, OCAL, OCCL, OHL, AVCAL, AVCDL
real(8) AVHL, AVPAL, AVCSL, AVCLL, LDENSMAGL, AVTL, LTHERMIFF, LVISCMAGL, LDTWRL
integer TransOnL, L1B, L2B, L3B, L4B, FFL, TEMPDXNL
real(8) TempDXL, OCABL, OHBL, AVCABL, AVHBL, AVPORBL, AVCSBL, AVCLBL, AVDENSBL, ALPHAMEIL
real(8) CRUSTMUL, MUSHSHEARL, CHECKMEIL, MAGVISCL, DIAML, LTHERMDIFFL, VISCMAGBL, HRTI0L
integer DMFNL, TRANSFL, tj, OUTCOUNTSM1L, maxmpi, topT, bottomT
real(8) DMFL, PRESSURERTIL, OCAT, OCCT, OHT, AVCAT, AVCDT, AVHT, AVPAT, AVCST, AVCLT
real(8) AVDENSILLT, AVTT, OVC1L, AVC1L, AVA1L, AVB1L, AVTSL, OC1T, AVC1T, AVB1T, AVA1T, AVTST
real(8) OC1, AvC1, Ava1, avb1, avts, oc1b, avc1b, ava1b, avb1b, avtsb
real(8) Oc1bl, ava1bl, avb1bl, avc1bl, avtsbl, maxmp

!isotopes
real(8) outisot, outci, Avia, oia, oial, oiccl, Avial, aviccl
real(8) oviccl, oiat, oicct, aviat, avicct
real(8) CompChoice

!isotopes - mantle
real(8) oimaL, oimcL, AvImaL,AvImcL, Oimat, Oimct
real(8) avimct, avimat,oima, oicc, oimc,avic, avimc
real(8) avima, outcim,outimt

real(8) DensFlagL, DensFlagM
integer denii, dmfnd

!new search
integer MT, MB, UT, UB
integer midD, upperD
integer oldj3, oldj4, oldL3, oldL4
integer moveupperD, movemidd
integer middlim, upperdlim, midDlim_LOWERONLY
real(8) diamU, diamM
real(8) rdmlm, rdmmu, rdmlu
integer ogmidlim,oguplim
real(8) lowD, middled
real(8) lowcrustm, midcrustm
real(8) lowRDM, midRDM
integer pori, depthli
integer porii, depthii

!Added for extraction efficiency
real(8) LEf, MEf
integer AmountLeaveL, AmountLeave
real(8) OcaTL, OCCTL, OHTL, OC1TL, OIATL, OIMATL, OICCTL, OIMCTL
real(8) AvCaTL, Avcdtl, avhtl, avc1tl, aviatl, avict, avimctl,avimatl
real(8) ava1tl, avb1tl, avtstl, avpatl, avcstl, avcltl, densmagtl, avttl
real(8)densmagLTL, avicctl
real(8) DensmagLT
real(8) dulc, tsm, cbm, a1m, b1m,c1m
integer dulcn

!added 12/06/23
real(8) runtime_erupt, time_erupt_s, buoy_fact
integer erupttimeF, top_shall, bottom_shall, outcountem1
real(8) ocae, occe, ohe, avcae, avcde, avhe, avpae, avcse, avcle, densmagle, avte
real(8) time_erupt_ka, oc1e, ava1e, avb1e, avc1e, avtse
character*5 FilecountEM1

!added 19/07/23
integer Melt_seg_flag, Melt_seg_flag_upper

!added 
integer EruptF

!added 17/01/25
real(8) low_DP, LDPL_OG, DPL_OG
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Set Variables

! Set some key variables to zero - if program fails, may need to set more of these
! if it is assumed they are initially zero

outcount=0
outcountSM=1
outcountsm1=1
OUTCOUNTSM1L=1
outcountem1=1
keff_dummy=0
TSill_dummy=0
velflag=1
runtime=0
timeoutcount=0
cumerror = 0
maxvel = 0
maxdiff = 0
outcomp=0
outcc=0
outci=0
outisot=0
outimt=0
outcim=0
TopFlagP=-50000000
BaseFlagP=500000000
AmountS=0
AmountI=0
AmountE=0
AILower=0
AELower=0
SN=0
oldL3=5000000
oldL4=5000000
oldj3=5000000
oldj4=5000000
EruptTimeF=0
runtime_erupt=0

ovtn=0
transn=0
checkMEI=0
transcount=0

do i = -n,n
	tb(i) = 0
	ta(i)=0
	hb(i) = 0
	ha(i)=0
	cb(i) = 0
	ca(i)=0
	cc(i)=0
	clb(i) = 0
	cla(i)=0
	csb(i) = 0
	csa(i) = 0
	pb(i) = 0
	pa(i) = 0
    ws(i) = 0
    wl(i) = 0
	icb(i) = 0
	ica(i) = 0
	icc(i) = 0
	icsb(i) = 0
	icsa(i) = 0
	iclb(i) = 0
	icla(i) = 0
	imb(i) = 0
	ima(i) = 0
	imsb(i) = 0
	imsa(i) = 0
	imlb(i) = 0
	imla(i) = 0
    imc(i) = 0
    pmag(i) = 0
	clostim(i) = -999
	advecttype(i) = 0
end do

! Read in from external file same folder as exe file

open(unit=20, file='A_DT_TE_Melt_Seg_NEW_1.txt')
read(20,*)low_DP
read(20,*)dummy !Main Choice - Flags!
read(20,*)dummy
read(20,*)MMorMFlag
read(20,*)dummy
read(20,*)EruptF
read(20,*)dummy
read(20,*)TransF
read(20,*)dummy
read(20,*)TransFL
read(20,*)dummy
read(20,*)WriteFL
read(20,*)dummy
read(20,*)maxpor
read(20,*)dummy
read(20,*)CompChoice
read(20,*)dummy
read(20,*)Melt_Seg_Flag
read(20,*)dummy
read(20,*)Melt_seg_flag_upper
read(20,*)dummy!transfer events!
read(20,*)dummy!main choices - mid crust transfer
read(20,*)dummy
read(20,*)time_erupt_ka
read(20,*)dummy
read(20,*)diamU
read(20,*)dummy
read(20,*)diamM
read(20,*)dummy
read(20,*)DMF
read(20,*)dummy
read(20,*)crustmu
read(20,*)dummy
read(20,*)Mef
read(20,*)dummy
read(20,*)Uef
read(20,*)dummy
read(20,*)DensFlagM
read(20,*)dummy
read(20,*)dummy!main choices - lower crust transfer
read(20,*)Tran_rate
read(20,*)dummy
read(20,*)TransWidth
read(20,*)dummy
read(20,*)diamL
read(20,*)dummy
read(20,*)DMFL
read(20,*)dummy
read(20,*)crustmuL
read(20,*)dummy
read(20,*)Lef
read(20,*)dummy
read(20,*)DensFlagL
read(20,*)dummy !background choices
read(20,*)dummy
read(20,*)buoy_fact
read(20,*)dummy
read(20,*)cpa
read(20,*)dummy
read(20,*)xabove
read(20,*)dummy
read(20,*)VolFlag
read(20,*)dummy
read(20,*)mpc
read(20,*)dummy
read(20,*)mc
read(20,*)dummy
read(20,*)Elsc
read(20,*)dummy
read(20,*)propL
read(20,*)dummy
read(20,*)JPBeta
read(20,*)dummy
read(20,*)tempDX
read(20,*)dummy
read(20,*)tempDXL
read(20,*)dummy
read(20,*)EmHMF
read(20,*)dummy !Physical properties!
read(20,*)dummy
read(20,*)a
read(20,*)dummy
read(20,*)b
read(20,*)dummy
read(20,*)alpha
read(20,*)dummy
read(20,*)beta
read(20,*)dummy
read(20,*)cee
!read(20,*)dummy
!read(20,*)drho
read(20,*)dummy
read(20,*)g
read(20,*)dummy
read(20,*)SLT
read(20,*)dummy
read(20,*)rho
read(20,*)dummy
read(20,*)min_melt_viscosity
read(20,*)dummy
read(20,*)max_melt_viscosity
!read(20,*)dummy
!read(20,*)bulk_viscosity
read(20,*)dummy
read(20,*)shear_viscosity
read(20,*)dummy
read(20,*)kt
read(20,*)dummy
read(20,*)Lf
read(20,*)dummy
read(20,*)Cp
read(20,*)dummy !Phase diagram 
read(20,*)dummy
read(20,*)Ae
read(20,*)dummy
read(20,*)dul
read(20,*)dummy
read(20,*)dulC
read(20,*)dummy
read(20,*)dummy ! top half 
read(20,*)TsU
read(20,*)dummy
read(20,*)cbu
read(20,*)dummy
read(20,*)A1u
read(20,*)dummy
read(20,*)B1u
read(20,*)dummy
read(20,*)C1u
read(20,*)dummy
read(20,*)dummy ! bottom half 
read(20,*)TsL
read(20,*)dummy
read(20,*)cbL
read(20,*)dummy
read(20,*)A1L
read(20,*)dummy
read(20,*)B1L
read(20,*)dummy
read(20,*)C1L
read(20,*)dummy ! 'below' depth of crust
read(20,*)dummy
read(20,*)TsM
read(20,*)dummy
read(20,*)cbM
read(20,*)dummy
read(20,*)A1M
read(20,*)dummy
read(20,*)B1M
read(20,*)dummy
read(20,*)C1M
read(20,*)dummy ! rest of crust
read(20,*)dummy
read(20,*)A2
read(20,*)dummy
read(20,*)B2
read(20,*)dummy
read(20,*)C2
read(20,*)dummy
read(20,*)aSiO2
read(20,*)dummy
read(20,*)bSiO2
read(20,*)dummy
read(20,*)cSiO2
read(20,*)dummy
read(20,*)dSiO2
read(20,*)dummy
read(20,*)DensSolA
read(20,*)dummy
read(20,*)DensSolB
read(20,*)dummy
read(20,*)DensLiqA
read(20,*)dummy
read(20,*)DensLiqB
read(20,*)dummy
read(20,*)densflag
read(20,*)dummy !Isotope properties
read(20,*)dummy
read(20,*)PartitionA
read(20,*)dummy
read(20,*)PartitionB
read(20,*)dummy
read(20,*)InitialI1
read(20,*)dummy
read(20,*)InitialI2
read(20,*)dummy
read(20,*)CloseTemp
read(20,*)dummy
read(20,*)StartTime
read(20,*)dummy !Sill properties
read(20,*)dummy 
read(20,*)sill_width
read(20,*)dummy
read(20,*)empexp
read(20,*)dummy
read(20,*)silliso1
read(20,*)dummy
read(20,*)silliso2
read(20,*)dummy
read(20,*)sill_comp
read(20,*)dummy
read(20,*)Emplacement_Depth1
read(20,*)dummy
read(20,*)Randdepth1
read(20,*)dummy
read(20,*)SillRandom1
read(20,*)dummy
read(20,*)Emplacement_Width1
read(20,*)dummy
read(20,*)Emplacement_Rate1
read(20,*)dummy
read(20,*)Emplacement_Pause1
read(20,*)dummy
read(20,*)Emplacement_Depth2
read(20,*)dummy
read(20,*)Randdepth2
read(20,*)dummy
read(20,*)SillRandom2
read(20,*)dummy
read(20,*)Emplacement_Width2
read(20,*)dummy
read(20,*)Emplacement_Rate2
read(20,*)dummy
read(20,*)Emplacement_Pause2
read(20,*)dummy
read(20,*)Emplacement_Depth3
read(20,*)dummy
read(20,*)Randdepth3
read(20,*)dummy
read(20,*)SillRandom3
read(20,*)dummy
read(20,*)Emplacement_Width3
read(20,*)dummy
read(20,*)Emplacement_Rate3
read(20,*)dummy
read(20,*)porodepth
read(20,*)dummy
read(20,*)accrete
read(20,*)dummy !General model properties
read(20,*)dummy
read(20,*)base
read(20,*)dummy
read(20,*)top
read(20,*)dummy
read(20,*)geotherm
read(20,*)dummy
read(20,*)Model_Time
read(20,*)dummy !Output data frequency and selected depth
read(20,*)dummy
read(20,*)Output_Interval
read(20,*)dummy
read(20,*)depthoutput
read(20,*)dummy
read(20,*)timeoutfreq
read(20,*)dummy
read(20,*)writetime
read(20,*)dummy !Solver/setup
read(20,*)dummy
read(20,*)keff_dummy
read(20,*)dummy
read(20,*)TSill_dummy
read(20,*)dummy
read(20,*)velflag
read(20,*)dummy
read(20,*)dzn
read(20,*)dummy
read(20,*)diffmin
read(20,*)dummy
read(20,*)diffexp
read(20,*)dummy
read(20,*)diffmax
read(20,*)dummy
read(20,*)tmultad1
read(20,*)dummy
read(20,*)tmultad2
read(20,*)dummy
read(20,*)tmultdiff
read(20,*)dummy !Mass conservation check
read(20,*)dummy
read(20,*)topflagU
read(20,*)dummy
read(20,*)baseflagU
read(20,*)dummy
read(20,*)topflagL
read(20,*)dummy
read(20,*)baseflagL
!read(20,*)dummy
!read(20,*)dtsill
close(unit=20)


!Open file to write everything from screen
open(unit=100, file='echo-screen.txt')


! Open txt files for various outputs 
open(unit=98, file='maxporosity-LOWER.txt')
open(unit=981, file='maxporosity-MID.txt')
open(unit=99, file='maxporosity-UPPER.txt')
open(unit=101, file='Total_Transferrable_Amount_Mid.txt')
open(unit=102, file='Total_Transferrable_Comp_Mid.txt')
open(unit=1012, file='Total_Transferrable_Iso_Mid.txt')
open(unit=1011, file='Total_Transferrable_Amount_Low.txt')
open(unit=1021, file='Total_Transferrable_Comp_Low.txt')
open(unit=1022, file='Total_Transferrable_Iso_Low.txt')
open(unit=103, file='A_MajorElements%.txt')
open(unit=1031, file='A_Isotopes_Conservation_crust%.txt')
open(unit=1032, file='A_Isotopes_Conservation_mantle%.txt')
open(unit=105, file='Total_Intruded_in_Shallow.txt')
open(unit=1051, file='Total_Intruded_Mid.txt')
open(unit=106, file='Total_Erupted.txt')
open(unit=1060, file='Total_Erupted_Comp.txt')
open(unit=107, file='LowerCheckN_Mid.txt')
open(unit=108, file='CompFlags.txt')
open(unit=109, file='Emplacement.txt')
open(unit=110, file='UpperCheckN_Shallow.txt')
open(unit=111, file='Magma_Buoy_Mid.txt')
open(unit=112, file='RTI_Buoy_Mid.txt')
open(unit=113, file='Tot_Buoy_Mid.txt')
open(unit=114, file='Crit_Buoy_Mid.txt')
open(unit=115, file='Magma_Buoy_Lower.txt')
open(unit=116, file='RTI_Buoy_Lower.txt')
open(unit=117, file='Tot_Buoy_Lower.txt')
open(unit=118, file='Crit_Buoy_Lower.txt')
open(unit=119, file='Magma_Buoy_Mid_C.txt')
open(unit=120, file='RTI_Buoy_Mid_C.txt')
open(unit=121, file='Tot_Buoy_Mid_C.txt')
open(unit=122, file='Crit_Buoy_Mid_C.txt')
open(unit=123, file='Magma_Buoy_Lower_C.txt')
open(unit=124, file='RTI_Buoy_Lower_C.txt')
open(unit=125, file='Tot_Buoy_Lower_C.txt')
open(unit=126, file='Crit_Buoy_Lower_C.txt')
open(unit=127, file='hdot_L.txt')
open(unit=128, file='hdot_M.txt')
open(unit=129, file='t1_L.txt')
open(unit=130, file='t1_M.txt')
open(unit=131, file='hmelt_L.txt')
open(unit=132, file='hmelt_M.txt')
open(unit=133, file='hbouy_L.txt')
open(unit=134, file='hbouy_M.txt')
open(unit=135, file='deltadens_L.txt')
open(unit=136, file='deltadens_M.txt')
open(unit=137, file='hRTI_L.txt')
open(unit=138, file='hRTI_M.txt')
open(unit=139, file='tauE_L.txt')
open(unit=140, file='tauE_M.txt')
open(unit=141, file='t2_L.txt')
open(unit=142, file='t2_M.txt')
open(unit=143, file='MagmaBuoyantVisc_L.txt')
open(unit=144, file='MagmaBuoyantVisc_M.txt')
open(unit=145, file='checkMEi_L.txt')
open(unit=146, file='checkMEi_M.txt')
open(unit=147, file='AmountLeave_vs_transfer.txt')
open(unit=148, file='AmountLeave_vs_transfer_lower.txt')
open(unit=213,file='checkDensL.txt')
open(unit=214,file='checkMu.txt')
open(unit=315,file='checkDensS.txt')
open(unit=316,file='checkSio2.txt')
open(unit=317,file='No_transfer_reason.txt')
open(unit=318,file='Depth_Intrusion_Mid.txt')
open(unit=319, file='A_bottom of upper crust.txt')
open(unit=320, file='A_bottom of mid crust.txt')
open(unit=3201, file='A_bottom of mid crust lower only.txt')
open(unit=321, file='A_oldj_.txt')
open(unit=322, file='A_j_.txt')
open(unit=323, file='A_oldl_.txt')
open(unit=324, file='A_l_.txt')
open(unit=325, file='A_diameteroflowT.txt')
open(unit=326, file='A_diameterofmidT.txt')
open(unit=327, file='A_lowTCV.txt')
open(unit=328, file='A_midTCV.txt')
open(unit=329, file='A_OGUPLIM.txt')
open(unit=330, file='A_OGMIDLIM.txt')

open(unit=390, file='Initial_temperature_mid.txt')
open(unit=391, file='Initial_tempertature_up.txt')




!! checks of the composition, viscosity and density.
do cacheck=0,1,0.05
    ndCheck=ndsio(cacheck,aSio2,bSio2,cSio2,dSio2)
    Sio2check=aSio2-bSiO2*tanh(-cSiO2*cacheck+dSiO2)
    densLCheck=((1-ndsio(cacheck,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(cacheck,aSio2,bSio2,cSio2,dSio2)*DensLiqB)
    densSCheck=((1-ndsio(cacheck,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(cacheck,aSio2,bSio2,cSio2,dSio2)*DensSolB)
    viscCheck=(10**(ndsio(cacheck,aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity)))*min_melt_viscosity
    write(213,*) caCheck,densLCheck
    write(214,*) caCheck,viscCheck
    write(315,*) cacheck,densSCheck
    write(316,*) cacheck,Sio2check
end do



! Initialising of more variables
drho = DensSolA - DensLiqB
Tl=C1U
Hs = Cp*TsU					!J
Hl = Cp*Tl + Lf				!J

hdot=0
hdotl=0
hLM1=0
hlm1l=0
TRANSON=0
transonl=0
trnsinprog=0

t1=0
t2=0
t1l=0
t2l=0
hrtiPL=0
hlm1L=0
hrtiP=0
hlm1=0
lambda=0
lambdaP=0
lambdaL=0
lambdaPL=0

BaseFlagP=BaseFlagL
TopFlagP=TopFlagL

DiamU=DiamU*1000
DiamM=DiamM*1000
DiamL=DiamL*1000


RDMLM=((DiamL/2)**2.0)/((DiamM/2)**2.0)

RDMMU=((DiamM/2)**2.0)/((DiamU/2)**2.0)

RDMLU=((DiamL/2)**2.0)/((DiamU/2)**2.0)

! Set sill intrusion temperature to be user defined or liquidus for chosen composition
! Tsill currently forced to be at liquidus as currently set SillPhi = 1
!if(TSill_dummy.gt.0)then
!    TSill = TSill_dummy
!else
TSill=Liquidus(sill_comp,Ae,A1U,B1U,C1U,A2,B2,C2)
BasDens = (ndsio(sill_comp,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+((1-ndsio(sill_comp,aSio2,bSio2,cSio2,dSio2))*DensLiqA)
!endif

write(*,*)'TSill = ',TSill-273.15,' deg C'
write(100,*)'TSill = ',TSill-273.15,' deg C'
write(*,*)'Sill density = ',BasDens,' kg/m3'
write(100,*)'Sill density = ',BasDens,' kg/m3'
write(*,*)'Reference density contrast = ',drho,' kg/m3'
write(100,*)'Reference density contrast = ',drho,' kg/m3'

Seed = 829376431

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Define Derived Parameters

perm = a*a*b

delta = sqrt((shear_viscosity+(4.0/3.0)*shear_viscosity)*perm/min_melt_viscosity)
tau = (1/(drho*g))*sqrt(min_melt_viscosity*(shear_viscosity+(4.0/3.0)*shear_viscosity)/perm)
omega = delta/tau

if(keff_dummy.gt.0)then
    keff = keff_dummy
else
    keff = kt*tau*(Tl-Tsu)/(rho*delta*delta*(cp*(Tl-Tsu) + Lf))
endif

stefan = Lf/(cp*(Tl-Tsu) + Lf)

Year = 365.25*24*60*60

write(*,*)'keff = ',keff
write(100,*)'keff = ',keff
write(*,*)'delta = ',delta,' m'
write(100,*)'delta = ',delta,' m'
write(*,*)'tau = ',tau/year,' yr'
write(100,*)'tau = ',tau/year,' yr'


! Sets frequency with which output files are written
OutRate = Output_Interval*1e3*Year/tau

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Define Modelling Paramters

dz = (Sill_Width*1000.0/delta)/dzn
DMFn = int((-DMF*1000)/(dz*delta))
tempDXn = int((tempDX*1000)/(dz*delta))
duln = int((-dul*1000)/(dz*delta))
dulCn = int((-dulC*1000)/(dz*delta))
tempDXnL = int((tempDXL*1000)/(dz*delta))
DMFnL = int((-DMFL*1000)/(dz*delta))
xaboven = int((xabove*1000)/(dz*delta))+1
midD=DMFnL
upperD=DMFn

write(319,*) 0.0, DMFn-xaboven
write(320,*) 0.0, DmFnL-xaboven


OGMIDLIM= DmFnL-xaboven

OGUPLIM= DMFn-xaboven





dtdiff = (dz**2)/(tmultdiff*keff)
dtadv = dz / tmultad1
!laxc = laxc*dz
!laxi = laxi*dz

!diffmin = diffmin*dz
!diffmax = diffmax*dz

if(diffmin.gt.diffmax)then 
    diffmin = diffmax
endif
!if(laxi.gt.dz)then 
!    laxi = dz
!endif

randdepthn1 = int((randdepth1*1000)/(dz*delta))
randdepthn2 = int((randdepth2*1000)/(dz*delta))
randdepthn3 = int((randdepth3*1000)/(dz*delta))

write(*,*)'dz = ',dz*delta,' m'
write(*,*)'dtdiff = ',dtdiff*tau/year,' yr'
write(*,*)'dtadv = ',dtadv*tau/year,' yr'
write(*,*)'dzn = ',dz,' -'
write(*,*)'dtndiff = ',dtdiff,' -'
write(*,*)'dtnadv = ',dtadv,' -'
write(100,*)'dz = ',dz*delta,' m'
write(100,*)'dtdiff = ',dtdiff*tau/year,' yr'
write(100,*)'dtadv = ',dtadv*tau/year,' yr'
write(100,*)'dzn = ',dz,' -'
write(100,*)'dtndiff = ',dtdiff,' -'
write(100,*)'dtnadv = ',dtadv,' -'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Define Non-Dimensionalised Parameters

z1 = int(-base*1000.0/(delta*dz))+1
z2 = int(-top*1000.0/(delta*dz))
d11 = z2 - int(Emplacement_Depth1*1000.0/(delta*dz))
d12 = z2 - int(Emplacement_Depth2*1000.0/(delta*dz))
d13 = z2 - int(Emplacement_Depth3*1000.0/(delta*dz))
d21 = d11
d22 = d12
d23 = d13

! Set initial emplacement depth 
emplacementi=d21
emplacementpi=z2
emplacementdi=z2
emplacementflagp=0
emplacementflagd=0	
sillstyle='initial depth'

time = (model_time*1e6*Year)/tau
ND_Geotherm = Geotherm*delta/(1000*(Tl-Tsu))

time_erupt_s=(time_erupt_ka*1e3*Year)/tau

!Dimensionless sill Paramters
SillNodes = int(Sill_Width*1000.0/(delta*dz))
SillTemp = (TSill - Tsu)/(Tl-Tsu)
SillPhi = 1.0
!Only works if Tsill = Tl otherwise will get Lf wrong and sill porosity wrong
SillEnth = (TSill*Cp+Lf - Hs)/(Hl-Hs)
SillComp = Sill_Comp

SillRate1 = (Emplacement_Rate1*1e6*Year)/tau
SillRate2 = (Emplacement_Rate2*1e6*Year)/tau
SillRate3 = (Emplacement_Rate3*1e6*Year)/tau

SillN1 = (Emplacement_Width1/Sill_Width)
SillN2 = (Emplacement_Width2/Sill_Width)
SillN3 = (Emplacement_Width3/Sill_Width)

SillCount1=0
SillCount2=0
SillCount3=0

weighttemp = weight

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Determine Starting Conditions And Find Node Corresponding to Output Depth
testflag=0
do i = z1,z2
	tb(i) = -i*dz*ND_Geotherm + (273.15-Tsu)/(Tl - Tsu)
    ta(i) = -i*dz*ND_Geotherm + (273.15-Tsu)/(Tl - Tsu)
! Only allocates enthalpy correctly if no melting so max T must be less than Ts
	hb(i) = enth_dim_to_nondim(Cp*temp_nondim_to_dim(tb(i),Tsu,Tl),Hs,Hl)
    
   
       Melt_seg(i)=1
    
        
    if (i.gt.duln) then 
        ca(i) = CbU
        cb(i) = CbU
        cc(i) = CbU
        clb(i) = 0
        csb(i) = CbU
        cla(i) = 0
        csa(i) = CbU
        Ts(i)=TsU
        if (Melt_seg_flag_UPPER.gt.0) then
            Melt_seg(i)=0
        end if
    else if ((i.le.duln).AND.(i.gt.dulCn)) then
        ca(i) = CbL
        cb(i) = CbL
        cc(i) = CbL
        clb(i) = 0
        csb(i) = CbL
        cla(i) = 0
        csa(i) = CbL
        Ts(i)=TsL
    else
        ca(i) = CbM
        cb(i) = CbM
        cc(i) = CbM
        clb(i) = 0
        csb(i) = CbM
        cla(i) = 0
        csa(i) = CbM
        Ts(i)=TsM
    end if
    
    pa(i) = 0
	pb(i) = 0
    ws(i) = 0
    wl(i) = 0
	icb(i) = InitialI1
    ica(i) = InitialI1
	icsb(i) = InitialI1
    icsa(i) = InitialI1
    icc(i) = InitialI1
	iclb(i) = 0
	imb(i) = InitialI2
    ima(i) = InitialI2
	imsb(i) = InitialI2
    imsa(i) = InitialI2
    !imla(i) = 0
	imlb(i) = 0
    pmag(i) = 0 
    imc(i) = InitialI2
    
	denscont(i) = 0
	density(i) = (ndsio(cb(i),aSio2,bSio2,cSio2,dSio2)*DensSolB)+((1-ndsio(cb(i),aSio2,bSio2,cSio2,dSio2))*DensSolA)
	
! Find node at which to write specific depth data
	if(((dz*i*delta/1000).gt.-depthoutput).AND.(testflag.lt.1)) then 
	    depthi = i
	    testflag=1
    endif
  
end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Output Initial Data

	open(1000,file='output0000.txt',status='old',err=20)
	close(1000,status='delete')
20	open(1000,file='output0000.txt',status='new')

	write(1000,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
	write(1000,'(22A20)') 'Distance (km)','Enthalpy (J)','Temperature (C)','Porosity','Melt Seg','P Mag','Solid Velocity (m/s)','Liquid Velocity (m/s)','Bulk Composition','Solid Composition','Liquid Composition',&
	&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
	&'Bulk Crust I','Solid Crust I','Liquid Crust I','Bulk Mantle I','Solid Mantle I','ts','Closure Time (My)','Density (kg/m3)'
	write(*,*), 'Initial Output Complete'
do i = z1,z2
	write(1000,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
	&F19.10,A1,F19.10,A1,F19.10,A14,F19.10)') &
	&dz*i*delta/1000,' ',hb(i),' ',(tb(i)*(Tl-Tsu)+Tsu)-273.15,' ',pb(i),' ',melt_seg(i),' ',pmag(i),' ',ws(i),' ',wl(i),' ',cb(i),' ',csb(i),' ',clb(i),&
	&' ',aSio2-bSiO2*tanh(-cSiO2*cb(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csb(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*clb(i)+dSiO2),&
	&' ',icb(i),' ',icsb(i),' ',iclb(i),' ',imb(i),' ',imsb(i),' ',ts(i),'        Undef ',density(i)
end do

	outcount = outcount + 1

	close(1000)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Bulk Solving Routine

!Open file to write all time-dependent data at maximum porosity
open(unit=50, file='output-time.txt')
write(50,'(19A20)') 'Time (My)','Depth (km)','Temperature (C)','Porosity','Bulk Composition','Solid Composition','Liquid Composition',&
&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
&'Bulk Crust I','Solid Crust I', 'Liquid Crust I', 'Bulk Mantle I','Solid Mantle I', 'Liquid Mantle I', 'Ratio Icb/Imb', 'Ratio Icl/Iml',&
&'Closure Time (My)'

!Open file to write time-dependent data at maximum porosity > SLT
open(unit=60, file='output-time-slt.txt')
write(60,'(19A20)') 'Time (My)','Depth (km)','Temperature (C)','Porosity','Bulk Composition','Solid Composition','Liquid Composition',&
&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
&'Bulk Crust I','Solid Crust I', 'Liquid Crust I', 'Bulk Mantle I','Solid Mantle I', 'Liquid Mantle I', 'Ratio Icb/Imb', 'Ratio Icl/Iml',&
&'Closure Time (My)'

!Open file to write time-dependent data at ALL porosity > SLT
open(unit=70, file='output-frequency-slt.txt')
write(70,'(19A20)') 'Time (My)','Depth (km)','Temperature (C)','Porosity','Bulk Composition','Solid Composition','Liquid Composition',&
&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
&'Bulk Crust I','Solid Crust I', 'Liquid Crust I', 'Bulk Mantle I','Solid Mantle I', 'Liquid Mantle I', 'Ratio Icb/Imb', 'Ratio Icl/Iml',&
&'Closure Time (My)'

!Open file to write time-dependent data at specific depth
open(unit=80, file='output-frequency-depth.txt')
write(80,'(18A20)') 'Time (My)','Temperature (C)','Porosity','Bulk Composition','Solid Composition','Liquid Composition',&
&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
&'Bulk Crust I','Solid Crust I', 'Liquid Crust I', 'Bulk Mantle I','Solid Mantle I', 'Liquid Mantle I', 'Ratio Icb/Imb', 'Ratio Icl/Iml',&
&'Closure Time (My)'

dsilltime = 0
freqcount=0

call pdvar(z1,z2,duln,dulcn,A1U,B1U,C1U,A1L,B1L,C1L,A1M,B1M,C1M,A1,B1,C1)
call S(z1,z2,duln,dulcn,TsU,TsL,TsM,Ts)

! Emplace first sill
write(*,*)'Sill emplacement depth controlled by ',sillstyle		
d2rand = d21	                   
call sillemplace(z1,z2,d11,d2rand,SillNodes,Accrete,&
&ha,ta,pa,ca,csa,cla,ica,icsa,icla,ima,imsa,imla,&
&hb,tb,pb,ws,wl,cb,cc,csb,clb,icb,icc,imc,icsb,iclb,imb,imsb,imlb,SillEnth,SillTemp,SillPhi,SillComp,SillIso1,&
&SillIso2,Ts,Tl,PartitionA,PartitionB,Ae,A1,B1,C1,A2,B2,C2,dz,delta,clostim,SN, A1U, B1U, C1U, TSU,runtime,tau,year,CloseTemp,melt_seg,pmag)	


dsilltime = 0
SillCount1 = SillCount1+1
TransOn=0
! no need for dt1 can replace with dt			     
dtl = dtdiff


! Start of runtime 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
do while (runtime.lt.time)
    
    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
    


seed=seed+runtime

! A check to make sure a sill is intruded below high melt fraction (EmHMF)
jEmH=5000000

do i=z1,z2
    if (pa(i).ge.EmHMF) then
        jEmH=i
        goto 99
    end if
end do




!Emplace Sill
99	if(SillCount1.lt.SillN1)then
		if (runtime.gt.real(SillCount1)*SillRate1-dtl) then
			if (runtime.lt.real(SillCount1)*SillRate1+dtl) then
	
    ! Find emplacement node
	            do i = z1,z2
    ! Based on porosity
                if((pa(i).gt.porodepth).AND.(emplacementflagp.lt.1))then
                    emplacementpi = i
                    emplacementflagp = 1
                endif
    ! Based on density
                if((BasDens.gt.density(i)).AND.(emplacementflagd.lt.1).AND.(densflag.gt.0))then
                    emplacementdi = i-1
                    emplacementflagd = 1
                endif
                
                  
                
                
                
    ! Choose deepest emplacement depth based on density or porosity, or use intial depth
                    if(abs(emplacementpi).gt.abs(emplacementdi))then
                    if(abs(emplacementpi).lt.1)then
                        emplacementi = d21
                        sillstyle='initial depth' 
                        multdepth = 1                  
                    else
                        emplacementi = emplacementpi
                        sillstyle='porosity'
                    endif
                else
                    if(abs(emplacementdi).lt.1)then
                        emplacementi = d21
                        sillstyle='initial depth' 
                        multdepth = 1                  
                    else
                        emplacementi = emplacementdi
                        sillstyle='density'
                    endif
                endif
	            end do
	            
                call RandomNumber(Random,randdepthn1,Seed,empexp)
!	            write(*,*)emplacementi,d21, emplacementdi,emplacementpi
	            
				if((porodepth.lt.1.000000001).OR.(densflag.gt.0))then
				     d2rand = emplacementi - random
				     				       
				else
                   if(accrete.gt.0)then
                     ! Overaccrete			    
				     d2rand = d21 - random
				     write(*,*)'Sill emplacement depth controlled by overaccretion around depth one'
				     write(100,*)'Sill emplacement depth controlled by overaccretion around depth one'
                 else
                     ! Underaccrete by setting emplacement depth just below previous sill plus some random difference
                     d2rand = d21 - SillCount1*SillNodes - random
                     write(*,*)'Sill emplacement depth controlled by underaccretion around depth one'
                     write(100,*)'Sill emplacement depth controlled by underaccretion around depth one'
                 endif 
                 endif  
    
    
                if (d2rand.gt.jEmH) then
                    d2rand=jEmh-1
                end if
                
                if (d2rand.ge.OGMIDLIM) then                     
                    d2rand=d21
                    sillstyle='No density - initial depth instead'
                end if
                
515               write(109,*) jEmH*dz*delta/1000!, d2rand*dz*delta/1000
                  write(*,*)'Sill emplacement depth controlled by ',sillstyle
                  write(100,*)'Sill emplacement depth controlled by ',sillstyle	
!              if(abs(emplacementpi).gt.abs(emplacementdi))then
!                  if(abs(emplacementpi).le.abs(d21))then
!                      emplacementi = d21
!                      sillstyle='initial depth' 
!                      multdepth = 1                  
!                  else
!                      emplacementi = emplacementpi
!                      sillstyle='porosity'
!                  endif
!              else
!                  emplacementi = emplacementdi
!                  sillstyle='density'
!              endif
!	            end do
!
!			    call RandomNumber(Random,randdepthn1,Seed,empexp)
!			    write(*,*)emplacementi,d21,emplacementdi,emplacementpi
!				if((porodepth.lt.1.000000001).OR.(densflag.gt.0))then
!				     d2rand = emplacementi - random
!				     write(*,*)'Sill emplacement depth controlled by ',sillstyle				       
!				else
!
!               if(accrete.gt.0)then
!                   ! Overaccrete			    
!				     d2rand = d21 - random
!				     write(*,*)'Sill emplacement depth controlled by overaccretion around depth one'
!                 else
!                     ! Underaccrete by setting emplacement depth just below previous sill plus some random difference
!                     d2rand = d21 - SillCount1*SillNodes - random
!                     write(*,*)'Sill emplacement depth controlled by underaccretion around depth one'
!                 endif 
!                 endif 
                   
	             call sillemplace(z1,z2,d11,d2rand,SillNodes,Accrete,&
		        &ha,ta,pa,ca,csa,cla,ica,icsa,icla,ima,imsa,imla,&
		        &hb,tb,pb,ws,wl,cb,cc,csb,clb,icb,icc,imc,icsb,iclb,imb,imsb,imlb,SillEnth,SillTemp,SillPhi,SillComp,SillIso1,&
			    &SillIso2,Ts,Tl,PartitionA,PartitionB,Ae,A1,B1,C1,A2,B2,C2,dz,delta,clostim,SN, A1U, B1U, C1U, TSU,runtime,tau,year,CloseTemp,melt_seg,pmag)
                 
                 ! Updating tracers of melt fraction
				
                 if ((oldj3.le.d2rand).AND.(oldj3.ne.5000000)) then 
                     oldj3=oldj3-SillNodes
                 end if
                 
                 if ((oldj4.le.d2rand).AND.(oldj4.ne.5000000)) then 
                     oldj4=oldj4-SillNodes
                 end if
                 
                 if ((oldL3.le.d2rand).AND.(oldL3.ne.5000000)) then 
                     oldL3=oldL3-SillNodes
                 end if
                 
                 if ((oldL4.le.d2rand).AND.(oldL4.ne.5000000)) then 
                     oldL4=oldL4-SillNodes
                 end if
                 
                 if (OGUPLIM.le.d2rand) then 
                    OGUPLIM=OGUPLIM-SillNodes
                 end if
                 
                 if (OGmidLIM.le.d2rand) then 
                    OGmidLIM=OGmidLIM-SillNodes
                 end if 
                 
                 
			     dsilltime = 0
			     SillCount1 = SillCount1+1
			     
			     ! Reset emplacement flag so finds new depth for next timestep
	             emplacementflagp=0
	             emplacementflagd=0	
                 emplacementpi=z2
                 emplacementdi=z2
	             
			end if
		end if
	end if

	if((SillCount1.ge.SillN1).AND.(SillCount2.lt.SillN2))then
		if (runtime.gt.real(SillCount1)*SillRate1+real(SillCount2)*SillRate2-dtl) then
			if (runtime.lt.real(SillCount1)*SillRate1+real(SillCount2)*SillRate2+dtl) then
				
			    call RandomNumber(Random,randdepthn2,Seed,empexp)
			    Seed = Seed + 2

!				if((porodepth.gt.0).AND.(emplacementi.lt.z2))then
!				   d2rand = emplacementi - random
!				
!				else
							    
                    if(accrete.gt.0)then
                ! Overaccrete			    
				        d2rand = d22 - random
				        write(*,*)'Sill emplacement depth controlled by overaccretion around depth two'
				        write(100,*)'Sill emplacement depth controlled by overaccretion around depth two'
                    else
                ! Underaccrete by setting emplacement depth just below previous sill plus some random difference
                        d2rand = d22 - random - (SillCount2*SillNodes)
                        write(*,*)'Sill emplacement depth controlled by underaccretion around depth two'
                        write(100,*)'Sill emplacement depth controlled by underaccretion around depth two'
                    endif 
!               endif
                if (d2rand.gt.jEmH) then
                    d2rand=jEmh-1
                end if
                
                if (d2rand.ge.OGMIDLIM) then 
                    d2rand=d21
                end if
                
516             write(109,*) jEmH*dz*delta/1000, d2rand*dz*delta/1000
                
                
	            call sillemplace(z1,z2,d12,d2rand,SillNodes,Accrete,&
					  &ha,ta,pa,ca,csa,cla,ica,icsa,icla,ima,imsa,imla,&
					  &hb,tb,pb,ws,wl,cb,cc,csb,clb,icb,icc,imc,icsb,iclb,imb,imsb,imlb,SillEnth,SillTemp,SillPhi,SillComp,SillIso1,&
					  &SillIso2,Ts,Tl,PartitionA,PartitionB,Ae,A1,B1,C1,A2,B2,C2,dz,delta,clostim,SN, A1U, B1U, C1U, TSU,runtime,tau,year,CloseTemp,melt_seg,pmag)
                
                ! Updating tracers of melt fraction
					  
                if ((oldj3.le.d2rand).AND.(oldj3.ne.5000000)) then 
                     oldj3=oldj3-SillNodes
                 end if
                 
                 if ((oldj4.le.d2rand).AND.(oldj4.ne.5000000)) then 
                     oldj4=oldj4-SillNodes
                 end if
                 
                 if ((oldL3.le.d2rand).AND.(oldL3.ne.5000000)) then 
                     oldL3=oldL3-SillNodes
                 end if
                 
                 if ((oldL4.le.d2rand).AND.(oldL4.ne.5000000)) then 
                     oldL4=oldL4-SillNodes
                 end if
                 
                 if (OGUPLIM.le.d2rand) then 
                    OGUPLIM=OGUPLIM-SillNodes
                 end if
                 
                 if (OGmidLIM.le.d2rand) then 
                    OGmidLIM=OGmidLIM-SillNodes
                 end if 
			    dsilltime = 0
			    SillCount2 = SillCount2+1
				
			end if
		end if
	end if
	
		if((SillCount1.ge.SillN1).AND.(SillCount2.ge.SillN2).AND.(SillCount3.lt.SillN3))then
		    if (runtime.gt.real(SillCount1)*SillRate1+real(SillCount2)*SillRate2+real(SillCount3)*SillRate3-dtl) then
			    if (runtime.lt.real(SillCount1)*SillRate1+real(SillCount2)*SillRate2+real(SillCount3)*SillRate3+dtl) then

			    call RandomNumber(Random,randdepthn3,Seed,empexp)
			    Seed = Seed + 2

!				if((porodepth.gt.0).AND.(emplacementi.lt.z2))then
!				   d2rand = emplacementi - random
				
!				else
				
                    if(accrete.gt.0)then
                ! Overaccrete			    
				        d2rand = d23 - random
				        write(*,*)'Sill emplacement depth controlled by overaccretion around depth three'
				        write(100,*)'Sill emplacement depth controlled by overaccretion around depth three'
                else
                ! Underaccrete by setting emplacement depth just below previous sill plus some random difference
                        d2rand = d23 - random - SillNodes - (SillCount3*SillNodes)
                        write(*,*)'Sill emplacement depth controlled by underaccretion around depth three'
                        write(100,*)'Sill emplacement depth controlled by underaccretion around depth three'
                    endif 
!                endif
                if (d2rand.gt.jEmH) then
                    d2rand=jEmh-1
                end if	
                
               if (d2rand.ge.OGMIDLIM) then                
                    d2rand=d21
                end if
                
517             write(109,*) jEmH*dz*delta/1000, d2rand*dz*delta/1000
                
	            call sillemplace(z1,z2,d13,d2rand,SillNodes,Accrete,&
					  &ha,ta,pa,ca,csa,cla,ica,icsa,icla,ima,imsa,imla,&
					  &hb,tb,pb,ws,wl,cb,cc,csb,clb,icb,icc,imc,icsb,iclb,imb,imsb,imlb,SillEnth,SillTemp,SillPhi,SillComp,SillIso1,&
					  &SillIso2,Ts,Tl,PartitionA,PartitionB,Ae,A1,B1,C1,A2,B2,C2,dz,delta,clostim,SN, A1U, B1U, C1U, TSU,runtime,tau,year,CloseTemp,melt_seg,pmag)
                
                ! Updating tracers of melt fraction
				if ((oldj3.le.d2rand).AND.(oldj3.ne.5000000)) then 
                     oldj3=oldj3-SillNodes
                 end if
                 
                 if ((oldj4.le.d2rand).AND.(oldj4.ne.5000000)) then 
                     oldj4=oldj4-SillNodes
                 end if
                 
                 if ((oldL3.le.d2rand).AND.(oldL3.ne.5000000)) then 
                     oldL3=oldL3-SillNodes
                 end if
                 
                 if ((oldL4.le.d2rand).AND.(oldL4.ne.5000000)) then 
                     oldL4=oldL4-SillNodes
                 end if	  				
                 
                 if (OGUPLIM.le.d2rand) then 
                    OGUPLIM=OGUPLIM-SillNodes
                 end if
                 
                 if (OGmidLIM.le.d2rand) then 
                    OGmidLIM=OGmidLIM-SillNodes
                 end if 
				SillCount3 = SillCount3+1
                dsilltime = 0
			end if
		end if
        end if

        
        
	!Update Enthalpy
	call enthsolve(z1,z2,dtl,dz,keff,stefan,ha,hb,tb,pb,ws,wl)
    
    
    !Find top and bottom of porosity.
    do i=z1,z2
        if ((pb(i).gt.0).OR.(ca(i).ne.cc(i)).OR.(enth_nondim_to_dim(ha(i),Hs,Hl).ge.(cp*ts(i))).OR.((ta(i)*(Tl-Tsu)+Tsu).ge.Ts(i))) then
            if (i.gt.TopFlagP) then
                TopFlagP=i
            end if
            if (i.lt.BaseFlagP) then
                BaseFlagP=i
                
            end if
        end if
    end do
    
   
       
    
	
	!Update Bulk Composition and Isotopes
	call compsolve(BaseFlagP, TopFlagP, z1,z2,dtl,dz,pa,pb,ca,cb,ws,wl,csb,clb,maxvel,maxdiff,cc,cerror,cerrorU,cerrorL,baseflagU,topflagU,baseflagL,topflagL,outcomp,outcc,&
                &diffexp,diffmax,diffmin,ica,icb,iclb,icsb,icc,ierror,outisot,outci,ima, imb,imlb,imsb,imc,imerror,outimt,outcim,CompChoice)
 
	!Update Porosity
	call porossolve(BaseFlagP,TopFlagP,ha,pa,ca,Hs,Hl,Lf,Cp,Tl,Ts,Ae,A1,B1,C1,A2,B2,C2,dum)
    
	!Update Temperature
	call tempsolve(z1,z2,ha,ta,pa,ca,Hs,Hl,Ts,Tsu,Tl,Lf,Cp)
	
	!Update Liquid and Solid Compositions
	do i = z1,z2
		if(pa(i).lt.0.00001) then
			cla(i) = 0.0
			csa(i) = ca(i)
			icla(i) = 0.0
			icsa(i) = ica(i)
			imla(i) = 0.0
			imsa(i) = ima(i)
		elseif(pa(i).gt.0.99999) then
			cla(i) = ca(i)
			csa(i) = 0.0
			icla(i) = ica(i)
			icsa(i) = 0.0
			imla(i) = ima(i)
			imsa(i) = 0.0
		else
			cla(i) = LiquidComp(temp_nondim_to_dim(ta(i),Tsu,Tl),ca(i),Ts(i),Ae,A1(i),B1(i),C1(i),A2,B2,C2)
			csa(i) = SolidComp(temp_nondim_to_dim(ta(i),Tsu,Tl),ca(i),Ae,A1(i),B1(i),C1(i),A2,B2,C2,Ts(i))
			icla(i) = IsoLiquid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
			icsa(i) = IsoSolid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)	
			imla(i) = IsoLiquid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
			imsa(i) = IsoSolid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
		end if
    end do
	
      !Update density
    
    do i=BaseflagP,TopflagP
        density(i) = pa(i)*((1-ndsio(cla(i),aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(cla(i),aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-pa(i))*((1-ndsio(csa(i),aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(csa(i),aSio2,bSio2,cSio2,dSio2)*DensSolB)
    end do

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! EVACUATIONS PORTION OF CODE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
! finding depth of shallow chamber 
    UT=5000000
    UB=5000000
    
    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
    
    
    do i=upperDlim,z2
        if ((pa(i).gt.0).AND.(pa(i+1).eq.0)) then
            UT=i
            goto 1031
        end if
    end do
    
1031 if (UT.eq.5000000) then 
        goto 1041
    end if
    
    do i=UT,z1,-1
        if ((pa(i).gt.0).AND.(pa(i-1).eq.0)) then
            UB=i
            UPPERD=i
            goto 1041
        end if
    end do
    
1041 if ((UT.eq.5000000).AND.(UB.eq.5000000)) then 
        UPPERD=DMFn
     end if
     
     
     ! this is checking that 2km beneath the deepest melt fracton in mid crust, there isn't another layer of melt   
     
    do i=upperD,z1,-1
        MOVEupperd=0
   
        do j=i-1,i-xaboven-1,-1
            if (pa(j).gt.0) then 
                MOVEupperd=1
                if (i.le.(z1+xaboven+1)) then
                    upperD=z1
                    midD=z1
                    goto 366
                end if
            end if
        end do
            
        if (MOVEupperd.eq.0) then
            upperD=i
            goto 367
        end if
            
        
    end do  
     
     
     
! finding the depth of the mid chamber 
367 MT=5000000
    MB=5000000
    
    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
    
    
    do i=midDlim, upperDlim-1
        if ((pa(i).gt.0).AND.(pa(i+1).eq.0)) then
            MT=i
            goto 1051
    end if
end do

    
1051 if (MT.eq.5000000) then
        goto 1061
    end if
   
    do i=MT,z1,-1
        if ((pa(i).gt.0).AND.(pa(i-1).eq.0)) then
            MB=i
            midD=i
            goto 1061
        end if
    end do
        
    
1061 if ((MT.eq.5000000).AND.(MB.eq.5000000)) then 
        midD=DMFnL
     end if
     
! this is checking that 2km beneath the deepest melt fracton in mid crust, there isn't another layer of melt   
     
    do i=midD,z1,-1
        MOVEMIDD=0
   
        do j=i-1,i-xaboven-1,-1
            if (pa(j).gt.0) then 
                MOVEMIDD=1
                if (i.le.(z1+xaboven+1)) then
                    midD=z1
                    goto 366
                end if
                
            end if
        end do
            
        if (MOVEMIDD.eq.0) then
            midD=i
            goto 366
        end if
            

    end do  
    

    
 ! find the melt fraction in the lower crust
    
   ! Counters needed
366 L1=5000000
    L2=5000000
    L3=5000000
    L4=5000000

    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
    
    
    midDlim_LOWERONLY = midDlim - 1
    
    if (TsL.ne.TsU) then
        do i=z1,midDlim-1
            if ((Ts(i).eq.TsU).AND.(Ts(i+1).eq.TsL)) then
                midDlim_LOWERONLY = i
                goto 1201
            end if
            
        end do
    end if
    
        
                
    
    
    ! find the top of the melt fraction - lower                        
            
                
1201    do i=midDlim_LOWERONLY,z1,-1
            if (pa(i).gt.0) then
                L3=i
                go to 205
            end if
        end do
    
        if (L3.eq.5000000) then
            pressureEL=0
            TotPressureL=0
            taueL=0
            LDPL=0
            
            t1L=0
            t2L=0
            hrtiPL=0
            hrtiL=0
            hlm1L=0
            lambdaL=0
            lambdaPL=0
            oldL3=5000000
            oldL4=5000000
            goto 202
        end if
        
    ! Find the bottom of the melt fraction - lower
    
205 do i=L3,z1,-1
        if ((pa(i).gt.0).AND.(pa(i-1).eq.0)) then
                L4=i
                
                goto 206
        end if
    end do
    
   206 if ((L3.eq.5000000).OR.(L4.eq.5000000)) then
            pressureEL=0
            TotPressureL=0
            taueL=0
            LDPL=0
            
            t1L=0
            t2L=0
            hrtiPL=0
            hrtiL=0
            hlm1L=0
            lambdaL=0
            lambdaPL=0
            oldL3=5000000
            oldL4=5000000
            goto 202
       end if   
       
     ! Find the top of the high melt fraction - lower
       
       do i=L3,L4,-1
           if ((pa(i).gt.cpa)) then
               L1=i
               goto 203
           end if
       end do
       
        if (L1.eq.5000000) then
            TRANSONL=1
            goto 204
            
        end if
    ! Find the bottom of high porosity
    
203 do i=L1,L4,-1
        if ((pa(i).ge.cpa).AND.(pa(i-1).lt.cpa)) then
                L2=i
                goto 204
        end if
    end do
        
204 if ((L1.eq.5000000).OR.(L2.eq.5000000)) then 
            TRANSONL=1
    end if
    
    
    
!! find out the overlying density of the crust. 
ovdenRL=0

do i=L3+1,L3+xaboven
    ovdenRL=ovdenRL+density(i)
end do

denRL=ovdenRL/((abs((L3+xaboven)-(L3+1))+1))
                
L3b=L3
                
do i=L3,L4,-1
if (density(i).ge.denRL) then
    ovdenRL=ovdenRL+density(i)
    denRL=ovdenRL/((abs((L3+xaboven)-(i))+1))
    L3b=i-1
    if (i.eq.L4) then
        pressureEL=0
        TotPressureL=0
        taueL=0
        LDPL=0
            
        t1L=0
        t2L=0
        hrtiPL=0
        hrtiL=0
        hlm1L=0
        lambdaL=0
        lambdaPL=0
        oldL3=5000000
        oldL4=5000000
        goto 202
    end if
else 
    go to 197
end if
end do

!! check that the high melt fraction layer is buoyant (if present). 

197 L2b=5000000
    L1B=5000000
    FFL=5000000

if (TRANSONL.eq.0) then
        do i=L1,L2,-1
            if (density(i).lt.denRL) then 
                L1B=i
               
                goto 225
            end if 
        end do 

        if (L1b.eq.5000000) then
            TRANSONL=1
    
        end if
end if

225 if (transonL.eq.0) then 

        do i=L1B,L2,-1 
                FFL=i
                if ((density(i).lt.denRL).AND.(density(i-1).ge.denRL)) then
                    L2B=i
                    if (L2b.gt.L1b) then
                         L2b=L1b
                    end if
                    goto 226
                end if
            end do
       226 if ((L2B.eq.5000000).AND.(FFL.lt.5000000).AND.(density(L2).lt.denRL)) then 
                L2B=L2
            end if
    
            if ((L2B.eq.5000000).OR.(L1B.eq.5000000)) then
                !! stop transfer
                TransonL=1
            end if 
end if 

! L1b is the top of the buoyant high melt fraction layer, L2b is the bottom of the buoyant high melt fraction layer. 
 
! L3b is the top of the buoyant magma layer, L4b is the bottom of the buoyant magma layer.

    L4b=L4
    do i=L3b,L4,-1
        if (density(i).ge.denRL) then 
            L4b=i+1
            if (L4b.gt.L3b) then
                L4b=L3b
            end if
            go to 198
        end if
    end do     
   
    !!Checking if this is new melt fraction or has previously been present
198    if ((oldL3.eq.5000000).OR.(oldL4.eq.5000000)) then
        oldL3=L3
        oldL4=L4
        goto 194
    else
        
        do i=L3,L4,-1
            do j=oldL3,oldL4,-1
                if (i.eq.j) then 
                    oldL3=L3
                    oldL4=L4
                    goto 194
                end if
            end do
        end do
    end if
    
    oldL3=L3
    oldL4=L4
    
    t1L=0
    t2L=0
    hrtiPL=0
    hlm1L=0
    lambdaL=0
    lambdaPL=0

    
!! Calculating the amount of melt (hlm) in metres
194 hmbL=(abs(L3b-L4b)+1)*dz*delta
   hlmL=(abs(L1b-L2b)+1)*dz*delta
!! Calculating the growth rate of hlm
    gRhlL=(hmbL-hlm1L)/(dtl*tau)
    !setting hl(t-1) for the next loop.
    hlm1L=hmbL
!! Caculating h dot
    hdotL=gRHlL
!! Calculating the amount of buoyant melt (hmb) ( in metres)
    
    if (L3.lt.OGMIDLIM) then 
        LOWd=DiamL
        LOWCrustM=crustmuL
    else if ((L3.ge.OGMIDLIM).AND.(L3.lt.OGUPLIM)) THEN 
        LOWd=DiamM
        LowCrustM=crustmu
    else
        LOWd=DiamU
        LowCrustM=crustmu
    end if 
    
if (WriteFL.gt.0) then    
    write(131,*) runtime*tau/(year*1000), hlmL, L1b,L2b
    write(127,*) runtime*tau/(year*1000), hdotL
    write(133,*) runtime*tau/(year*1000), hmbL,L3b,L4b
end if 
 
  !!!!!!!! average things for transfer  - only happens if there is high melt fraction 
    if (TRANSONL.eq.0) then
            !1) find average comp and enth
                OcaL=0
                OccL=0
                OhL=0
                OvC1L=0
                OiaL=0
                oiccL=0
                OimaL=0
                oimcL=0
                
                do i=L2b,L1b
                    OcaL=OcaL+ca(i)
                    OhL=OhL+ha(i)
                    OccL=OccL+cc(i)
                    OvC1L=OvC1L+C1(i)
                    OiaL=OiaL+ica(i)
                    OiccL=OiccL+icc(i)
                    OimaL=OimaL+ima(i)
                    OimcL=OimcL+imc(i)
                end do 
                
                AvcaL=OcaL/(abs(L2b-L1b)+1)
                AvcdL=OccL/(abs(L2b-L1b)+1)
                AvhL=OhL/(abs(L2b-L1b)+1)
                AvIaL=OiaL/(abs(L2b-L1b)+1)
                AviccL=OiccL/(abs(L2b-L1b)+1)
                AvImaL=OimaL/(abs(L2b-L1b)+1)
                AvImcL=OimcL/(abs(L2b-L1b)+1)
                
                AvC1L=OvC1L/(abs(L2b-L1b)+1)
                
                if (abs(AvC1L-C1U).lt.abs(AvC1L-C1L)) then 
                    AvA1L=A1u
                    AvB1L=B1u
                    AvC1L=C1u
                    AvtsL=TsU
                else
                    AvA1L=A1l
                    AvB1L=B1l
                    AvC1L=C1l
                    AvtsL=tsl
                end if
                
                    
                
            !2) find porosity
               call porosshift(AvhL,AvpaL,AvcaL,Hs,Hl,Lf,Cp,Tl,AvTsL,Ae,AvA1L,AvB1L,AvC1L,A2,B2,C2)  
               
            !3) find composition of liquid and solid
               AvcsL=0
               AvclL=(1.0/AvpaL)*(AvcaL)
               if (AvclL.gt.1) then
                   AvclL=1.0
               end if
               
            !4) find average density of magmaL
               LdensMagL=AvpaL*((1-ndsio(AvclL,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(AvclL,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-AvpaL)*((1-ndsio(AvcsL,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(AvcsL,aSio2,bSio2,cSio2,dSio2)*DensSolB)      
               AvTL=AvA1L*AvclL**2+AvB1L*AvclL+AvC1L
              
               
               AvTL=temp_dim_to_nondim(AvTL,Tsu,Tl)       
    
    
    
    
    
!!!!!!!! finding the critical overpressure (DPL)

   
   
   LthermDiffL=kt/(Cp*LdensMagL)
   
   LViscMagL = (10**(ndsio(AvClL,aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity)))*min_melt_viscosity
   

   
   
   LdTWRL= (abs((temp_nondim_to_dim(ta(L1b),Tsu,Tl)-273.15)-(temp_nondim_to_dim(ta(L1b+tempDxnL),Tsu,Tl)-273.15)))&
            &/((abs(L1b-(L1b+tempDxnL))+1)*dz*delta)
   
   LDPL=((2.0*sqrt(3.0))/(JPbeta))**(2.0/5.0)*(sqrt((LthermDiffL*LViscMagL)/(pi_8))*((cp*LdTWRL)/(lf))*Elsc**2)**(2.0/5.0)
   
       if (LDPL<low_DP) then 
            LDPL_OG = LDPL
            LDPL = low_DP
       else
           LDPL_OG = LDPL
        end if
    
   
   !! outputs critical overpresure for the melt to leave the magma reservoir
        if (WriteFL.gt.0) then
            !! write pressure here
            write(118,*) runtime*tau/(year*1000), LDPL_OG
        end if   
   
    end if
 !!end of evaluating high melt fraction
    
    
    
!!!!!! finding the pressure from the buoyant melt layer alone (pressureE)

OCABL=0
OhbL=0
Oc1bl=0

!find average comp and enth - buoyant layer
do i=L3b,L4b,-1     
    OcabL=ocabL+ca(i)
    OhbL=OhbL+ha(i)
    oc1bl=Oc1bl+c1(i)
end do
AvcabL=ocabL/(abs(L3b-L4b)+1)
AvhbL=OhbL/(abs(L3b-L4b)+1)
Avc1bl=oc1bl/(abs(L3b-L4b)+1)

if ((abs(Avc1bl-C1U)).lt.(abs(Avc1bl-C1L))) then
    AvA1bL=A1U
    AvB1bl=B1U
    AvC1bl=C1U
    AvTsbl=TsU
else
    AvA1bL=A1L
    AvB1bl=B1L
    AvC1bl=C1L
    AvTsbl=TsL
end if



!find porosity - buoyant layer 
call porosshift(AvhbL,AvporbL,AvcabL,Hs,Hl,Lf,Cp,Tl,AvTsbL,Ae,AvA1bL,AvB1bL,AvC1bL,A2,B2,C2) 

!find comp - buoyant layer
AvcsbL=0
AvclbL=(1/AvporbL)*AvcabL

if (avclbL.gt.1.0) then
    AvclbL=1.0
end if


!find density - buoyant layer
AvdensbL=AvporbL*((1-ndsio(AvclbL,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(AvclbL,aSio2,bSio2,cSio2,dSio2)*DensLiqB)&
&+(1-AvporbL)*((1-ndsio(AvcsbL,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(AvcsbL,aSio2,bSio2,cSio2,dSio2)*DensSolB)   


if (L3b.eq.L4b) then
    AvdensBL=density(L3b)
    AvporBL=pa(L3b)
    AvclbL=cla(L3b)
end if 

if (AvdensBL.gt.denRL) then
    pressureEL=0
    PressureRTIL=0
    TotPressureL=0
    taueL=0
            
    t1L=0
    t2L=0
    hrtiPL=0
    hlm1L=0
    lambdaL=0
    lambdaPL=0
    oldL3=5000000
    oldL4=5000000
    goto 202
end if

!work out excess pressure
pressureEL=BUOY_FACT*(denRL-AvdensBL)*g*(abs(L3b-L4b)+1)*dz*delta

!write out excess pressure
if (WriteFL.gt.0) then
    !! write pressure here 
    write(115,*) runtime*tau/(year*1000), pressureEL
end if     

!Calculate viscosity of buoyant magma layer - if it is less than CPA mush shear viscosity from keller 
!This is only used a prediction of lambda to calculate the unconfined time, which is no longer used because
! of the findings of Booth et al 2024, as unconfined time << total time to evacuation. 
alphaMEIL=-((1.0/cpa)*log(max_melt_viscosity/LowCrustM))

MushShearL=LowCrustM*exp(-alphaMeiL*AvporBL)
ViscMagBL = (10**(ndsio(AvClbL,aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity)))*min_melt_viscosity

if (checkMEIL.eq.0) then
    checkMEIL=1
    write(145,*) alphaMEIL
    
    do i=1,11
        write(145,*) LowCrustM*exp(-alphaMeiL*(i-1)/10)
    end do
end if



if (AvPorBL.lt.cpa) then
    MagViscL=MushShearL
else
    MagViscL=max_melt_viscosity
end if




!!!!! calculating lambda
!!! cancel this out if needed
if (lambdaL.lt.LowD) then
    if (t1L.eq.0) then
        lambdaL=lambdaPL
        else 
        lambdaL=lambdaPL+((4.0*PI_8)/(2.88))*hdotL*dtl*tau*((LowCrustM/MagViscL)**(1.0/3.0))
    end if   
    !t2=0
!end if 

!!!!! lambda<D - which means chance of brittle failure
!if (lambda.lt.LowD) then 
    !update t1
    t1L=t1L+dtl*tau
    !!brittle failure check!!      
    !previous lambda
    lambdaPL=lambdaL
end if
!!!!
    !if (TransF.gt.0) then 
      !  if (pressureE.lt.DPL) then
     !       goto 202
    !    else
   !         TransK=1
            !failure due to brittle.
  !          goto 78
 !       end if
  !  end if 
    
!!! lamba=D
!else if (lambda.ge.LowD) then 
    !write(318,*) t1, lambda


! Calculate hRTI growth    
    tauEL=(6.0*pi_8*LowCrustM)/(LowD*g*(denRL-AvdensBL))
    
    
    
    hRTI0L=PROPL*LowD
    
    if (t2L.eq.0) then
        hRTIL=hRTI0L*exp(t2L/tauEL)
    else 
        
    
        hRTIL=hrtiPL*exp((dtl*tau)/tauEL)
    end if
    
    if (hRTIL.ge.(abs(L3b-midD)*delta*dz)) then
        hRTIL=abs(L3b-midD)*delta*dz
    end if 
    
    hRTIPL=hRTIL

   
    
    t2L=t2L+dtl*tau
    
!!! brittle + rti    
    pressureRTIL=BUOY_FACT*hRTIL*(denRL-AvdensBL)*g
 
    TotPressureL=pressureRTIL+pressureEL
    if (WriteFL.gt.0) then
        !WRITE PRESSURE HERE
        write(116,*) runtime*tau/(year*1000), pressureRTIL
        write(117,*) runtime*tau/(year*1000), TotPressureL
    end if 
        if (TransFL.gt.0) then 
            if ((pressureRTIL+pressureEL).lt.LDPL) then
                goto 202
            end if 
        end if


!end if



!Calculations of evacuations as criteria for evacuation is met

178 if ((TRANSONL.gt.0).OR.(TransFL.eq.0)) then 
         goto 202
    end if

     ! !!!! HERE - TRANSNODES CHANGES TO MATCH UP WITH SUBROUTINE
          
    
    AmountLeaveL=int((abs(L1b-L2b))*LEf)
    
    write(148,*) runtime*tau/(year*1000), abs(L1b-L2b), AmountLeaveL
    
    !find averages transferred
                !1) find average comp and enth
                OcaTL=0
                OccTL=0
                OhTL=0
                Oc1TL=0
                OiaTL=0
                OimaTL=0
                OiccTL=0
                OimcTL=0
                oimaTL=0
                
                
            
                do i=L1b-AmountLeaveL,L1b
                    OcaTL=OcaTL+ca(i)
                    OhTL=OhTL+ha(i)
                    OccTL=OccTL+cc(i)
                    Oc1TL=Oc1TL+c1(i)
                    oiaTL=oiaTL+ica(i)
                    oiccTL=oiccTL+icc(i)
                    oimcTL=oimcTL+imc(i)
                    oimaTL=oimaTL+ima(i)
                end do 
                
                AvcaTL=OcaTL/(AmountLeaveL+1)
                AvcdTL=OccTL/(AmountLeaveL+1)
                AvhTL=OhTL/(AmountLeaveL+1)
                Avc1TL=Oc1TL/(AmountLeaveL+1)
                aviaTL=oiaTL/(AmountLeaveL+1)
                AviccTL=oiccTL/(AmountLeaveL+1)
                AvimcTL=oimcTL/(AmountLeaveL+1)
                AvimaTL=oimaTL/(AmountLeaveL+1)
                
                if ((abs(avc1TL-c1u)).lt.(abs(avc1TL-c1l))) then
                        AvA1TL=A1U
                        AvB1TL=B1U
                        AvC1TL=C1U
                        AvtsTL=TsU
                else
                      AvA1TL=A1L
                      AvB1TL=B1L
                      AvC1TL=C1L
                      AvtsTL=TsL
                end if
    
            
                        
                
            !2) find porosity
               call porosshift(AvhTL,AvpaTL,AvcaTL,Hs,Hl,Lf,Cp,Tl,AvTsTL,Ae,AvA1TL,AvB1TL,AvC1TL,A2,B2,C2)  
               
            !3) find composition of liquid and solid
               AvcsTL=0
               AvclTL=(1.0/AvpaTL)*(AvcaTL)
               if (AvclTL.gt.1) then
                   AvclTL=1.0
               end if
               
            !4) find average density of magma
               densMagLTL=AvpaTL*((1-ndsio(AvclTL,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(AvclTL,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-AvpaTL)*((1-ndsio(AvcsTL,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(AvcsTL,aSio2,bSio2,cSio2,dSio2)*DensSolB)      
               AvTTL=AvA1TL*AvclTL**2+AvB1TL*AvclTL+AvC1TL
              
               
               AvTTL=temp_dim_to_nondim(AvTTL,Tsu,Tl)       
            
    ! Calculate depth that magma is evacuated to in the mid-crust
          if (DensFLagL.gt.0) then 
          
                deni=500000
                pori=500000
            
                do i=L3+xaboven+1,z2
            
                    if (DensMagLTL.ge.density(i)) then
                        deni=i
                        goto 800
                    end if
                end do
        
800             do i=L3+xaboven+1,z2
                    
                    if (pa(i).gt.porodepth) then 
                        pori=i
                        goto 803
                    end if
                    
                end do
                
803             if (deni.le.pori) then 
                    DepthLi=deni
                    write(100,*) runtime*tau/(year*1000), 'Lower transfer depth chosen by density'
                else
                    DepthLi=pori
                    write(100,*) runtime*tau/(year*1000), 'Lower transfer depth chosen by porosity'
                end if
                
                
                
                if (DensMagLTL.ge.DenRL) then 
                    write(317,*) runtime*tau/(year*1000), 'not buoyant enough'
                    write(100,*) runtime*tau/(year*1000), 'not buoyant enough - lower'
                    goto 202
                else if ((depthli.eq.500000).or.(depthli.eq.z2).OR.(depthli.ge.OGUPLIM)) then 
                    Tdepth=DMFnL
                else
                    Tdepth=depthli
                end if
                
          else
              Tdepth=DMFNL
          end if
          
               
          if (LowD.eq.diamL) then
                if (Tdepth.lt.OGMIDLIM) then 
                    lowRDM=1.0
                else if ((Tdepth.ge.OGMIDLIM).AND.(Tdepth.lt.OGUPLIM)) THEN 
                    lowRDM=RDMLM
                else
                    lowRDM=RDMLU
                end if 
          else if (LowD.eq.diamM) then 
                if (Tdepth.lt.OGMIDLIM) then 
                    lowRDM=0
                else if ((Tdepth.ge.OGMIDLIM).AND.(Tdepth.lt.OGUPLIM)) THEN 
                    lowRDM=1.0
                else
                    lowRDM=RDMMU
                end if
          else
              if (Tdepth.lt.OGMIDLIM) then 
                    lowRDM=0
                else if ((Tdepth.ge.OGMIDLIM).AND.(Tdepth.lt.OGUPLIM)) THEN 
                    lowRDM=0
                else
                    lowRDM=1.0
                end if
          end if
          
              
               
                
                
    
            !outputs  
                write(1011,*) runtime*tau/(year*1000),  (AmountLeaveL+1)*dz*delta/1000, ((AmountLeaveL+1)*dz*delta/1000)*((((LowD/2)/1000)**2)*PI_8) 
                write(1021,*) runtime*tau/(year*1000), L1b*dz*delta/1000,  aSio2-bSiO2*tanh(-cSiO2*AvcaTL+dSiO2) 
                write(1022,*) runtime*tau/(year*1000), AvIaTL,AvImaTL
                write(318,*) runtime*tau/(year*1000), TDepth*dz*delta/1000
                write(390,*) runtime*tau/(year*1000), (ta(TDepth)*(Tl-Tsu)+Tsu)-273.15
                          
                    if(outcountSM1L.lt.10) then
		            write(filecountSM1L,'(A1,A1,A1,A1,I1)') '0','0','0','0',outcountSM1L
		            else if (outcountSM1L.lt.100) then
		            write(filecountSM1L,'(A1,A1,A1,I2)') '0','0','0',outcountSM1L
		            else if (outcountSM1L.lt.1000) then
		            write(filecountSM1L,'(A1,A1,I3)') '0','0',outcountSM1L
                    else if (outcountSM1L.lt.10000) then
		            write(filecountSM1L,'(A1,I4)') '0',outcountSM1L
		            else
		            write(filecountSM1L,'(I5)') outcountSM1L
		            end if
		            outcountSM1L = outcountSM1L + 1

                  
                    
         
                    open(415,file='LowerCrust_Failure_Values_1'//filecountSM1L//'.txt',status='old',err=191)
	                close(415,status='delete')
191                 open(415,file='LowerCrust_Failure_Values_1'//filecountSM1L//'.txt',status='new')      
                    write(415,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(415,*) 'Total pressure (Pa)'
                    write(415,*) TotPressureL
                    write(415,*) 'Overlying density (kg/m^3)'
                    write(415,*) DenRL
                    write(415,*) 'Crtical Pressure (Pa)'
                    write(415,*) LDPL
                    write(415,*) 'h - km'
                    write(415,*) hmbL/1000
                    write(415,*) 'hRTI -km'
                    write(415,*) hRTIL/1000
                    write(415,*) 't1 (ka)'
                    write(415,*) t1L/(year*1000)
                    write(415,*) 't2 (ka)'
                    write(415,*) t2L/(year*1000)
                    write(415,*) 'RDM'
                    write(415,*) lowRDM
                    write(415,*) 'Tdepth'
                    write(415,*) Tdepth
                    write(415,*) 'OGmid'
                    write(415,*) OGMIDLIM
                    write(415,*) 'OGup'
                    write(415,*) OGupLIM
                    write(415,*) 'crustmu'
                    write(415,*) Lowcrustm
                    close(415)
                
                
        
            
                    open(416,file='LowerCrust_Failure_Values_2'//filecountSM1L//'.txt',status='old',err=181)
	                close(416,status='delete')
181                 open(416,file='LowerCrust_Failure_Values_2'//filecountSM1L//'.txt',status='new')      
                    write(416,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(416,*) 'Density of ALL the HMF magma layer'
                    write(416,*) LDensMagL
                    write(416,*) 'Density of evacuated HMF magma layer'
                    write(416,*) DensMagLTL
                    write(416,*) 'h density'
                    write(416,*) AvdensbL
                    write(416,*) 'Thermal Diffusivity - mag'
                    write(416,*) LthermDiffL
                    write(416,*) 'Viscosity -mag'
                    write(416,*) LViscMagL
                    write(416,*) 'Temperature Gradient'
                    write(416,*) LdtWRL
                    write(416,*) 'lambda'
                    write(416,*) lambdaL
                    write(416,*) 'hdot'
                    write(416,*) hdotL
                    write(416,*) 'Average Porosity - Buoyant'
                    write(416,*) AvporBL
                    write(416,*) 'MagVisc'
                    write(416,*) MagViscL
                    close(416)
  
                
                !! Transferred overview
                    open(418,file='temp-lower'//filecountSM1L//'.txt',status='old',err=121)
	                close(418,status='delete')
121                  open(418,file='temp-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(418,*) (ta(i)*(Tl-Tsu)+Tsu)-273.15
                    end do
                    close(418)
                    
                    open(419,file='porosity-lower'//filecountSM1L//'.txt',status='old',err=131)
	                close(419,status='delete')
131                  open(419,file='porosity-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(419,*) pa(i)
                    end do
                    close(419)                    
                    
                    open(420,file='comp-lower'//filecountSM1L//'.txt',status='old',err=141)
	                close(420,status='delete')
141                  open(420,file='comp-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(420,*) ca(i)
                    end do
                    close(420) 
                    
                    open(421,file='density-lower'//filecountSM1L//'.txt',status='old',err=151)
	                close(421,status='delete')
151                  open(421,file='density-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(421,*) density(i)
                    end do
                    close(421)                      
                    
                    open(422,file='isotopes-crust-lower'//filecountSM1L//'.txt',status='old',err=152)
	                close(422,status='delete')
152                 open(422,file='isotopes-crust-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(422,*) ica(i)
                    end do
                    close(422) 
                    
                    open(424,file='isotopes-mantle-lower'//filecountSM1L//'.txt',status='old',err=154)
	                close(424,status='delete')
154                 open(424,file='isotopes-mantle-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(424,*) ima(i)
                    end do
                    close(424)   
                    
                    open(425,file='isotopes-ratio-lower'//filecountSM1L//'.txt',status='old',err=155)
	                close(425,status='delete')
155                 open(425,file='isotopes-ratio-lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(425,*) ica(i)/ima(i)
                    end do
                    close(425)     
                    
                    
                    open(423,file='closure times - lower'//filecountSM1L//'.txt',status='old',err=153)
	                close(423,status='delete')
153                   open(423,file='closure times - lower'//filecountSM1L//'.txt',status='new') 
                    do i=L2b,L1b
                        write(423,*) clostim(i), runtime*tau/year/1E6-clostim(i)
                    end do
                    close(423)    
                    
                    
                    
!! Transferred averages        
                    open(417,file='Average_of_HMF_layer_lower'//filecountSM1L//'.txt',status='old',err=111)
	                close(417,status='delete')
111                  open(417,file='Average_of_HMF_layer_lower'//filecountSM1L//'.txt',status='new')      
                    write(417,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(417,*) 'Average density'
                    write(417,*) LDensMagL
                    write(417,*) 'Average comp'
                    write(417,*) AvcaL
                    write(417,*) 'Average crust isotopes'
                    write(417,*) AvIaL
                    write(417,*) 'Average mantle isotopes'
                    write(417,*) AvImaL
                    write(417,*) 'Average temp'
                    write(417,*) (AvTL*(Tl-Tsu)+Tsu)-273.15
                    write(417,*) 'Average porosity'
                    write(417,*) AvpaL
                    write(417,*) 'Average cl'
                    write(417,*) AvclL
                    close(417)
            
                    open(430,file='Averages_transferred_lower'//filecountSM1L//'.txt',status='old',err=112)
	                close(430,status='delete')
112                 open(430,file='Averages_transferred_lower'//filecountSM1L//'.txt',status='new')      
                    write(430,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(430,*) 'Average density'
                    write(430,*) DensMagLTL
                    write(430,*) 'Average comp'
                    write(430,*) AvcaTL
                    write(430,*) 'Average crust isotopes'
                    write(430,*) AvIaTL
                    write(430,*) 'Average mantle isotopes'
                    write(430,*) AvImaTL
                    write(430,*) 'Average temp'
                    write(430,*) (AvTTL*(Tl-Tsu)+Tsu)-273.15
                    write(430,*) 'Average porosity'
                    write(430,*) AvpaTL
                    write(430,*) 'Average cl'
                    write(430,*) AvclTL
                    close(430)
                            

                
                write(1051,*) runtime*tau/(year*1000), (AmountLeaveL+1)*lowRDM*dz*delta/1000, ((AmountLeaveL+1)*dz*delta/1000)*((((LowD/2)/1000)**2)*PI_8)
                
                ! subroutine for evacuation
                call lowershiftmass(L1b,L1b-AmountLeaveL,z1,z2,oldj3,oldj4,ovtn,AmountLeaveL+1,lowRDM,Tdepth,ha,ta,&
                    &pa,ca,cc,csa,cla,ica,icc,icsa,icla,ima,imc,&
                    &imsa,imla,hb,tb,pb,cb,csb,clb,icb,icsb,iclb,imb,imsb,imlb,clostim,density,&
                    & PartitionA,PartitionB,AvcaTL,AvcdTL,AvTTL,AvhTL,AvclTL,AvcsTL,AvpaTL,DensMagLTL,&
                    &runtime,year,tau,baseflagP,topflagP,dz,delta,asio2,bsio2,csio2,dsio2, &
                    &AvA1TL, AvB1TL, AvC1TL, AvTsTL,A1,B1,C1,Ts,AvIaTL,AvImaTL,AvIccTL,AvImcTL,Tl,Tsu,Closetemp,midD,OGMIDLIM,OGUPLIM,d21,melt_seg,pmag)
               

                !Reset values
               ! pressureEL=0
                !TotPressureL=0
                PressureRTIL=0
                taueL=0
                LDPL=5e20
                TRANSONL=0
                t1L=0
                t2L=0
                hrtiPL=0
                hlm1L=0
                lambdaL=0
                lambdaPL=0    
                oldL3=5000000
                oldL4=5000000
    
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ! Evacuations to the shallow crust/eruption
202 TransonL=0
    j1=5000000
    j2=5000000
    j3=5000000
    j4=5000000
    
    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
 
! find the top of melt fraction 
    do i=upperDlim-1, midDlim, -1
        if (pa(i).gt.0) then
            j3=i
            goto 105
        end if
    end do 
    
    if (j3.eq.5000000) then
            pressureE=0
            TotPressure=0
            taue=0
            DPL=0
            
            t1=0
            t2=0
            hrtiP=0
            hrti=0
            hlm1=0
            lambda=0
            lambdaP=0
            oldj3=5000000
            oldj4=5000000
            goto 102
    end if
   
! find the bottom of melt fraction
105  do i=j3,midDlim,-1
        if ((pa(i).gt.0).AND.(pa(i-1).eq.0)) then
                j4=i
                goto 106
        end if
    end do
        
106 if ((j3.eq.5000000).OR.(j4.eq.5000000)) then
            pressureE=0
            TotPressure=0
            taue=0
            DPL=0
            
            t1=0
            t2=0
            hrtiP=0
            hrti=0
            hlm1=0
            lambda=0
            lambdaP=0
            oldj3=5000000
            oldj4=5000000
            goto 102
    end if   
    
  
 ! find the top of the high melt fraction layer    
  do i=j3,j4,-1
        
        if ((pa(i).ge.cpa))then
            j1=i
            goto 103
        end if
        
        
        
    end do
    
        if (j1.eq.5000000) then
            TRANSON=1
            goto 104
            
        end if  
    
    
    ! Find the bottom of high melt fraction layer
    
103 do i=j1,j4,-1
        if ((pa(i).ge.cpa).AND.(pa(i-1).lt.cpa)) then
                j2=i
                goto 104
        end if
    end do
        
104 if ((j1.eq.5000000).OR.(j2.eq.5000000)) then 
            TRANSON=1

            
    end if    
    
    
    !! find out the overlying density of the crust. 
    ovdenR=0

    do i=j3+1,j3+xaboven
        ovdenR=ovdenR+density(i)
    end do

    denR=ovdenR/((abs((j3+xaboven)-(j3+1))+1))
                
    j3b=j3
                
    do i=j3,j4,-1
    if (density(i).ge.denR) then
        ovdenR=ovdenR+density(i)
        denR=ovdenR/((abs((j3+xaboven)-(i))+1))
        j3b=i-1
        if (i.eq.j4) then
            pressureE=0
            TotPressure=0
            taue=0
            DPL=0
            
            t1=0
            t2=0
            hrtiP=0
            hrti=0
            hlm1=0
            lambda=0
            lambdaP=0
            oldj3=5000000
            oldj4=5000000
            goto 102
        end if
    else 
        go to 97
    end if
    end do
    

!! check that the high melt fraction layer is buoyant (if present). 

97 j2b=5000000
j1B=5000000
FF=5000000

if (TRANSON.eq.0) then
        do i=j1,j2,-1
            if (density(i).lt.denR) then 
                j1B=i
                goto 125
            end if 
        end do 

        if (j1b.eq.5000000) then
            TRANSON=1
    
        end if
end if

125 if (transon.eq.0) then 

        do i=j1B,j2,-1 !!! here change
                FF=i
                if ((density(i).lt.denR).AND.(density(i-1).ge.denR)) then
                    j2B=i
                    if (j2b.gt.j1b) then
                         j2b=j1b
                    end if
                    goto 126
                end if
            end do
        126 if ((j2B.eq.5000000).AND.(FF.lt.5000000).AND.(density(j2).lt.denR)) then 
                j2B=j2
            end if
    
            if ((j2B.eq.5000000).OR.(j1B.eq.5000000)) then
                !! stop transfer
                Transon=1
            end if 
end if 

! j1b is the top of the buoyant high melt fraction layer, j2b is the bottom of the buoyant high melt fraction layer. 
 
! j3b is the top of the buoyant magma layer, j4b is the bottom of the buoyant magma layer.

    j4b=j4
    do i=j3b,j4,-1
        if (density(i).ge.denR) then 
            j4b=i+1
            if (j4b.gt.j3b) then
                j4b=j3b
            end if
            go to 98
        end if
    end do 

 ! check if melt fraction layer is new or the same layer as previous time step
98  if ((oldj3.eq.5000000).OR.(oldj4.eq.5000000)) then
        oldj3=j3
        oldj4=j4
        goto 94
    else
        
        do i=j3,j4,-1
            do j=oldj3,oldj4,-1
                if (i.eq.j) then 
                    oldj3=j3
                    oldj4=j4
                    goto 94
                end if
            end do
        end do
    end if
    
    oldj3=j3
    oldj4=j4

    t1=0
    t2=0
    hrtiP=0
    hlm1=0
    lambda=0
    lambdaP=0
    
!! Calculating the amount of melt (hlm) in metres
94 hmb=(abs(j3b-j4b)+1)*dz*delta
   hlm=(abs(j1b-j2b)+1)*dz*delta
!! Calculating the growth rate of hlm
    gRhl=(hmb-hlm1)/(dtl*tau)
    !setting hl(t-1) for the next loop.
    hlm1=hmb
!! Caculating h dot
    hdot=gRHl
!! Calculating the amount of buoyant melt (hmb) ( in metres)
    
    if (j3.lt.OGMIDLIM) then 
        MIDdled=DiamL
        MidCrustM=crustmuL
    else if ((j3.ge.OGMIDLIM).AND.(j3.lt.OGUPLIM)) THEN 
        MIDdled=DiamM
        MidCrustM=crustmu
    else
        MIDdled=DiamU
        MidCrustM=crustmu
    end if
 
 if (WriteFL.gt.0) then
    write(132,*) runtime*tau/(year*1000), hlm, j1b,j2b
    write(128,*) runtime*tau/(year*1000), hdot
    write(134,*) runtime*tau/(year*1000), hmb,j3b,j4b
 end if
 
 
  !!!!!!!! average things for transfer  - only happens if there is high melt fraction 
    if (TRANSON.eq.0) then
            !1) find average comp and enth
                Oca=0
                Occ=0
                Oh=0
                Oc1=0
                Oia=0
                Oima=0
                Oicc=0
                Oimc=0
                oima=0
                
                do i=j2b,j1b
                    Oca=Oca+ca(i)
                    Oh=Oh+ha(i)
                    Occ=Occ+cc(i)
                    Oc1=Oc1+c1(i)
                    oia=oia+ica(i)
                    oicc=oicc+icc(i)
                    oimc=oimc+imc(i)
                    oima=oima+ima(i)
                end do 
                
                Avca=Oca/(abs(j2b-j1b)+1)
                Avcd=Occ/(abs(j2b-j1b)+1)
                Avh=Oh/(abs(j2b-j1b)+1)
                Avc1=Oc1/(abs(j2b-j1b)+1)
                avia=oia/(abs(j2b-j1b)+1)
                Avic=oicc/(abs(j2b-j1b)+1)
                Avimc=oimc/(abs(j2b-j1b)+1)
                Avima=oima/(abs(j2b-j1b)+1)
                
                if ((abs(avc1-c1u)).lt.(abs(avc1-c1l))) then
                        AvA1=A1U
                        AvB1=B1U
                        AvC1=C1U
                        Avts=TsU
                else
                      AvA1=A1L
                      AvB1=B1L
                      AvC1=C1L
                      Avts=TsL
                end if
    
            
                        
                
            !2) find porosity
               call porosshift(Avh,Avpa,Avca,Hs,Hl,Lf,Cp,Tl,AvTs,Ae,AvA1,AvB1,AvC1,A2,B2,C2)  
               
            !3) find composition of liquid and solid
               Avcs=0
               Avcl=(1.0/Avpa)*(Avca)
               if (Avcl.gt.1) then
                   Avcl=1.0
               end if
               
            !4) find average density of magma
               densMagL=Avpa*((1-ndsio(Avcl,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(Avcl,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-Avpa)*((1-ndsio(Avcs,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(Avcs,aSio2,bSio2,cSio2,dSio2)*DensSolB)      
               AvT=AvA1*Avcl**2+AvB1*Avcl+AvC1
              
               
               AvT=temp_dim_to_nondim(AvT,Tsu,Tl)       
    
    
    
    
    
!!!!!!!! finding the critical overpressure (DPL)

   
   
   thermDiffL=kt/(Cp*densMagL)
   
   ViscMagL = (10**(ndsio(AvCl,aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity)))*min_melt_viscosity
   

   
   
   dTWRL= (abs((temp_nondim_to_dim(ta(j1b),Tsu,Tl)-273.15)-(temp_nondim_to_dim(ta(j1b+tempDxn),Tsu,Tl)-273.15)))&
            &/((abs(j1b-(j1b+tempDxn))+1)*dz*delta)
   
   DPL=((2.0*sqrt(3.0))/(JPbeta))**(2.0/5.0)*(sqrt((thermDiffL*ViscMagL)/(pi_8))*((cp*dTWRL)/(lf))*Elsc**2)**(2.0/5.0)
   
   if (DPL<low_DP) then
        DPL_OG = DPL
        DPL = low_DP
   else
       DPL_OG = DPL
    
   end if
    
   
   !! outputs critical overpresure for the melt to leave the magma reservoir
        if (WriteFL.gt.0) then
            !! write pressure here
            write(114,*) runtime*tau/(year*1000), DPL_OG
        end if   
   
    end if
 !!end of evaluating high melt fraction
    
    
    
!!!!!! finding the pressure from the buoyant melt layer alone (pressureE)

OCAB=0
Ohb=0
OC1b=0


!find average comp and enth - buoyant layer
do i=j3b,j4b,-1     
    Ocab=ocab+ca(i)
    Ohb=Ohb+ha(i)
    oc1b=oc1b+c1(i)
end do
Avcab=ocab/(abs(j3b-j4b)+1)
Avhb=Ohb/(abs(j3b-j4b)+1)
Avc1b=oc1b/(abs(j3b-j4b)+1)

if ((abs(Avc1b-c1U)).lt.(abs(avc1b-c1l))) then 
    AvA1b=A1U
    AvB1b=B1U
    AvC1b=C1U
    AvTsb=Tsu
else
    AvA1b=A1L
    AvB1b=B1L
    AvC1b=C1L
    AvTsb=Tsu
end if



!find porosity - buoyant layer 
call porosshift(Avhb,Avporb,Avcab,Hs,Hl,Lf,Cp,Tl,AvTsb,Ae,AvA1b,AvB1b,AvC1b,A2,B2,C2) 

!find comp - buoyant layer
Avcsb=0
Avclb=(1/Avporb)*Avcab

if (avclb.gt.1.0) then
    Avclb=1.0
end if


!find density - buoyant layer
Avdensb=Avporb*((1-ndsio(Avclb,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(Avclb,aSio2,bSio2,cSio2,dSio2)*DensLiqB)&
&+(1-Avporb)*((1-ndsio(Avcsb,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(Avcsb,aSio2,bSio2,cSio2,dSio2)*DensSolB)   


if (j3b.eq.j4b) then
    AvdensB=density(j3b)
    AvporB=pa(j3b)
    Avclb=cla(j3b)
end if 

if (AvdensB.gt.denR) then
    pressureE=0
    PressureRTI=0
    TotPressure=0
    taue=0
            
    t1=0
    t2=0
    hrtiP=0
    hlm1=0
    lambda=0
    lambdaP=0
    oldj3=5000000
    oldj4=5000000
    goto 102
end if
!work out excess pressure
pressureE=BUOY_FACT*(denR-AvdensB)*g*(abs(j3b-j4b)+1)*dz*delta

!write out excess pressure
if (WriteFL.gt.0) then
    write(111,*) runtime*tau/(year*1000), pressureE
end if     

!Calculate viscosity of buoyant magma layer - if it is less than CPA mush shear viscosity from keller 
!This is only used a prediction of lambda to calculate the unconfined time, which is no longer used because
! of the findings of Booth et al 2024, as unconfined time << total time to evacuation. 
alphaMEI=-((1.0/cpa)*log(max_melt_viscosity/MidCrustM))

MushShear=MidCrustM*exp(-alphaMei*AvporB)
ViscMagB = (10**(ndsio(AvClb,aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity)))*min_melt_viscosity

if (checkMEI.eq.0) then
    checkMEI=1
    write(146,*) alphaMEI
    
    do i=1,11
        write(146,*) MidCrustM*exp(-alphaMei*(i-1)/10)
    end do
end if



if (AvPorB.lt.cpa) then
    MagVisc=MushShear
else
    MagVisc=max_melt_viscosity
end if



!!!!! calculating lambda
!!! cancel this out if needed
if (lambda.lt.MIDdleD) then
    if (t1.eq.0) then
        lambda=lambdaP
        else 
        lambda=lambdaP+((4.0*PI_8)/(2.88))*hdot*dtl*tau*((MidCrustM/MagVisc)**(1.0/3.0))
    end if   
    !t2=0
!end if 

!!!!! lambda<D - which means chance of brittle failure
!if (lambda.lt.MIDdleD) then 
    !update t1
    t1=t1+dtl*tau
    !!brittle failure check!!      
    !previous lambda
    lambdaP=lambda
end if
!!!!
    !if (TransF.gt.0) then 
      !  if (pressureE.lt.DPL) then
     !       goto 102
    !    else
   !         TransK=1
            !failure due to brittle.
  !          goto 78
 !       end if
  !  end if 
    
!!! lamba=D
!else if (lambda.ge.MIDD) then 
    !write(318,*) t1, lambda



! Calculating hrti growth    
    tauE=(6.0*pi_8*MidCrustM)/(MIDdleD*g*(denR-AvdensB))
    
    
    
    hRTI0=PROPL*MIDdleD
    
    if (t2.eq.0) then
        hRTI=hRTI0*exp(t2/tauE)
    else 
        
    
        hRTI=hrtiP*exp((dtl*tau)/tauE)
    end if
    
    
    if (hRTI.ge.(abs(j3b-upperD)*delta*dz)) then
        hRTI=abs(j3b-upperD)*delta*dz
    end if 
    
    hRTIP=hRTI

   
    
    t2=t2+dtl*tau
    
!!! brittle + rti    
    pressureRTI=BUOY_FACT*hRTI*(denR-AvdensB)*g
    TotPressure=pressureRTI+pressureE
    if (WriteFL.gt.0) then
        !! WRITE PRESSURE HERE
        write(112,*) runtime*tau/(year*1000), pressureRTI
        write(113,*) runtime*tau/(year*1000), TotPressure
    end if 
        if (TransF.gt.0) then 
            if ((pressureRTI+pressureE).lt.DPL) then
                goto 102
            else 
                TransK=2
                !! failure due to brittle+RTI
                goto 78    
            end if 
        end if


!end if




!!!!!!!!!!!!!!!evacuation occurs!!!!!!!!!!!!!!
    
78 if (TRANSON.gt.0) then 
         goto 102
    end if

    
    
    
    !! If excess pressure is gt or equal to DPL - let it continue, if not end.
      
  if (TransF.gt.0) then
      
            AmountLeave=int((abs(j1b-j2b))*MEf)
             !!WRITE CODE AMOUNTLEAVE AND ABSJ1B
            write(147,*)  runtime*tau/(year*1000),abs(j1b-j2b),AmountLeave
            
            
            !1) find average comp and enth
                OcaT=0
                OccT=0
                OhT=0
                Oc1T=0
                OiaT=0
                OimaT=0
                OiccT=0
                OimcT=0
                oimaT=0
                
                
            
                do i=j1b-AmountLeave,j1b
                    OcaT=OcaT+ca(i)
                    OhT=OhT+ha(i)
                    OccT=OccT+cc(i)
                    Oc1T=Oc1T+c1(i)
                    oiaT=oiaT+ica(i)
                    oiccT=oiccT+icc(i)
                    oimcT=oimcT+imc(i)
                    oimaT=oimaT+ima(i)
                end do 
                
                AvcaT=OcaT/(AmountLeave+1)
                AvcdT=OccT/(AmountLeave+1)
                AvhT=OhT/(AmountLeave+1)
                Avc1T=Oc1T/(AmountLeave+1)
                aviaT=oiaT/(AmountLeave+1)
                AvicT=oiccT/(AmountLeave+1)
                AvimcT=oimcT/(AmountLeave+1)
                AvimaT=oimaT/(AmountLeave+1)
                
                if ((abs(avc1T-c1u)).lt.(abs(avc1T-c1l))) then
                        AvA1T=A1U
                        AvB1T=B1U
                        AvC1T=C1U
                        AvtsT=TsU
                else
                      AvA1T=A1L
                      AvB1T=B1L
                      AvC1T=C1L
                      AvtsT=TsL
                end if
    
            
                        
                
            !2) find porosity
               call porosshift(AvhT,AvpaT,AvcaT,Hs,Hl,Lf,Cp,Tl,AvTsT,Ae,AvA1T,AvB1T,AvC1T,A2,B2,C2)  
               
            !3) find composition of liquid and solid
               AvcsT=0
               AvclT=(1.0/AvpaT)*(AvcaT)
               if (AvclT.gt.1) then
                   AvclT=1.0
               end if
               
            !4) find average density of magma
               densMagLT=AvpaT*((1-ndsio(AvclT,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(AvclT,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-AvpaT)*((1-ndsio(AvcsT,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(AvcsT,aSio2,bSio2,cSio2,dSio2)*DensSolB)      
               AvTT=AvA1T*AvclT**2+AvB1T*AvclT+AvC1T
              
               
               AvTT=temp_dim_to_nondim(AvTT,Tsu,Tl)       
            
            
            
            
            
            
            
            
            
           
                
      ! Find depth of intrusion into the shallow crust
      
      
            if (DensFLagM.gt.0) then 
          
                denii=500000
                porii=500000
            
                do i=j3+xaboven+1,z2
            
                    if (DensMagLT.ge.density(i)) then
                        denii=i
                        goto 801
                    end if
                end do
        
801             do i=j3+xaboven+1,z2
                    if (pa(i).gt.porodepth) then 
                        porii=i
                        goto 802
                    end if
                end do
                
802             if (denii.le.porii) then 
                    Depthii=denii
                    write(100,*) runtime*tau/(year*1000),'Transfer depth to shallow crust choosen by density'
                else
                    Depthii=porii
                    write(100,*) runtime*tau/(year*1000),'Transfer depth to shallow crust choosen by porosity'
                end if
                
                if (DensMagLT.ge.DenR) then 
                    write(317,*) runtime*tau/(year*1000), 'not buoyant enough-upper'
                    write(100,*) runtime*tau/(year*1000),'no upper transfer'
                    goto 102
                else if ((Depthii.eq.500000).or.(Depthii.eq.z2)) then 
                    DMFnd=DMFn
                else
                    DMFnd=Depthii
                end if
                
                    
                
          else
              DMFNd=DMFN
          end if
       
        if (MIDdleD.eq.diamL) then
                if (DMFnd.lt.OGMIDLIM) then 
                    MIDRDM=1.0
                else if ((DMFnd.ge.OGMIDLIM).AND.(DMFnd.lt.OGUPLIM)) THEN 
                    MIDRDM=RDMLM
                else
                    MIDRDM=RDMLU
                end if 
        else if (MIDdleD.eq.diamM) then 
                if (DMFnd.lt.OGMIDLIM) then 
                    MIDRDM=0
                else if ((DMFnd.ge.OGMIDLIM).AND.(DMFnd.lt.OGUPLIM)) THEN 
                    MIDRDM=1.0
                else
                    MIDRDM=RDMMU
                end if
        else
                if (DMFnd.lt.OGMIDLIM) then 
                    MIDRDM=0
                else if ((DMFnd.ge.OGMIDLIM).AND.(DMFnd.lt.OGUPLIM)) THEN 
                    MIDRDM=0
                else
                    MIDRDM=1.0
                end if
         end if

            
            
              !outputs              
                    if(outcountSM1.lt.10) then
		            write(filecountSM1,'(A1,A1,A1,A1,I1)') '0','0','0','0',outcountSM1
		            else if (outcountSM1.lt.100) then
		            write(filecountSM1,'(A1,A1,A1,I2)') '0','0','0',outcountSM1
		            else if (outcountSM1.lt.1000) then
		            write(filecountSM1,'(A1,A1,I3)') '0','0',outcountSM1
                    else if (outcountSM1.lt.10000) then
		            write(filecountSM1,'(A1,I4)') '0',outcountSM1
		            else
		            write(filecountSM1,'(I5)') outcountSM1
		            end if
		            outcountSM1 = outcountSM1 + 1
                    
            if (TransK.eq.1) then    
                    
                    open(215,file='Brittle_val'//filecountSM1//'.txt',status='old',err=17)
	                close(215,status='delete')
17                  open(215,file='Brittle_val'//filecountSM1//'.txt',status='new')      
                    write(215,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(215,*) 'pe'
                    write(215,*) PressureE
                    write(215,*) 'Overlying density'
                    write(215,*) DenR
                    write(215,*) 'Crtical Pressure'
                    write(215,*) DPL
                    write(215,*) 'h'
                    write(215,*) hmb
                    close(215)
                    
            else if (TransK.eq.2) then
                    open(215,file='BrittleRTI_val'//filecountSM1//'.txt',status='old',err=19)
	                close(215,status='delete')
19                  open(215,file='BrittleRTI_val'//filecountSM1//'.txt',status='new')      
                    write(215,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(215,*) 'Total pressure'
                    write(215,*) TotPressure
                    write(215,*) 'Overlying density'
                    write(215,*) DenR
                    write(215,*) 'Crtical Pressure'
                    write(215,*) DPL
                    write(215,*) 'h - magma - km'
                    write(215,*) hmb/1000
                    write(215,*) 'hRTI -km'
                    write(215,*) hRTI/1000
                    write(215,*) 't1'
                    write(215,*) t1/(year*1000)
                    write(215,*) 't2'
                    write(215,*) t2/(year*1000)
                    write(215,*) 'RDM'
                    write(215,*) MIDRDM
                    write(215,*) 'DMFND'
                    write(215,*) DMFND
                    write(215,*) 'OGmid'
                    write(215,*) OGMIDLIM
                    write(215,*) 'OGup'
                    write(215,*) OGupLIM
                    write(215,*) 'crustmu'
                    write(215,*) midcrustm
                    write(215,*) 'middleD'
                    write(215,*) middleD
                    close(215)
                    
                
                
            end if
            
                    open(216,file='DPL_val'//filecountSM1//'.txt',status='old',err=18)
	                close(216,status='delete')
18                  open(216,file='DPL_val'//filecountSM1//'.txt',status='new')      
                    write(216,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(216,*) 'Density of ALL the HMF layer'
                    write(216,*) DensMagL
                    write(216,*) 'Density of evacuated HMF layer'
                    write(216,*) DensMagLT
                    write(216,*) 'buoyant magma density'
                    write(216,*) Avdensb
                    write(216,*) 'Thermal Diffusivity'
                    write(216,*) thermDiffL
                    write(216,*) 'Viscosity'
                    write(216,*) ViscMagL
                    write(216,*) 'Temperature Gradient'
                    write(216,*) dtWRL
                    write(216,*) 'Crtical Pressure'
                    write(216,*) DPL
                    write(216,*) 'lambda'
                    write(216,*) lambda
                    write(216,*) 'hdot'
                    write(216,*) hdot
                    write(216,*) 'Average Porosity - Buoyant'
                    write(216,*) AvporB
                    write(216,*) 'MagVisc'
                    write(216,*) MagVisc
                    close(216)

!! Transferred overview
                    open(218,file='temp'//filecountSM1//'.txt',status='old',err=12)
	                close(218,status='delete')
12                  open(218,file='temp'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(218,*) (ta(i)*(Tl-Tsu)+Tsu)-273.15
                    end do
                    close(218)
                    
                    open(219,file='porosity'//filecountSM1//'.txt',status='old',err=13)
	                close(219,status='delete')
13                  open(219,file='porosity'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(219,*) pa(i)
                    end do
                    close(219)                    
                    
                    open(220,file='comp'//filecountSM1//'.txt',status='old',err=14)
	                close(220,status='delete')
14                  open(220,file='comp'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(220,*) ca(i)
                    end do
                    close(220) 
                    
                    open(221,file='density'//filecountSM1//'.txt',status='old',err=15)
	                close(221,status='delete')
15                  open(221,file='density'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(221,*) density(i)
                    end do
                    close(221)    
                    
                   open(222,file='isotopes-crust'//filecountSM1//'.txt',status='old',err=16)
	                close(222,status='delete')
16                  open(222,file='isotopes-crust'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(222,*) ica(i)
                    end do
                    close(222)
                    
                    open(224,file='isotopes-mantle'//filecountSM1//'.txt',status='old',err=29)
	                close(224,status='delete')
29                  open(224,file='isotopes-mantle'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(224,*) ima(i)
                    end do
                    close(224)
                    
                    open(227,file='isotopes-ratio'//filecountSM1//'.txt',status='old',err=30)
	                close(227,status='delete')
30                  open(227,file='isotopes-ratio'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(227,*) ica(i)/ima(i)
                    end do
                    close(227)
                    
                    open(223,file='closure times'//filecountSM1//'.txt',status='old',err=185)
	                close(223,status='delete')
185                  open(223,file='closure times'//filecountSM1//'.txt',status='new') 
                    do i=j2b,j1b
                        write(223,*) clostim(i), runtime*tau/year/1E6-clostim(i)
                    end do
                    close(223)
                    
                    
!! Transferred averages        
                    open(217,file='Average_of_HMF_layer'//filecountSM1//'.txt',status='old',err=11)
	                close(217,status='delete')
11                  open(217,file='Average_of_HMF_layer'//filecountSM1//'.txt',status='new')      
                    write(217,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(217,*) 'Average density'
                    write(217,*) DensMagL
                    write(217,*) 'Average comp'
                    write(217,*) Avca
                    write(217,*) 'Average Crust Isotopes'
                    write(217,*) AvIa
                    write(217,*) 'Average Mantle Isotopes'
                    write(217,*) AvIma
                    write(217,*) 'Average temp'
                    write(217,*) (AvT*(Tl-Tsu)+Tsu)-273.15
                    write(217,*) 'Average porosity'
                    write(217,*) Avpa
                    write(217,*) 'Average cl'
                    write(217,*) Avcl
                    close(217)
                    
                    open(228,file='Averages_transferred'//filecountSM1//'.txt',status='old',err=110)
	                close(228,status='delete')
110                  open(228,file='Averages_transferred'//filecountSM1//'.txt',status='new')      
                    write(228,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(228,*) 'Average density'
                    write(228,*) DensMagLT
                    write(228,*) 'Average comp'
                    write(228,*) AvcaT
                    write(228,*) 'Average Crust Isotopes'
                    write(228,*) AvIaT
                    write(228,*) 'Average Mantle Isotopes'
                    write(228,*) AvImaT
                    write(228,*) 'Average temp'
                    write(228,*) (AvTT*(Tl-Tsu)+Tsu)-273.15
                    write(228,*) 'Average porosity'
                    write(228,*) AvpaT
                    write(228,*) 'Average cl'
                    write(228,*) AvclT
                    close(228)
                    
                    
                    
                    
            !! Transferred amount and composition
              write(101,*) runtime*tau/(year*1000), (AmountLeave+1)*dz*delta/1000, ((AmountLeave+1)*dz*delta/1000)*((((MIDdleD/2)/1000)**2)*PI_8) 
              write(102,*) runtime*tau/(year*1000), j1b*dz*delta/1000, aSio2-bSiO2*tanh(-cSiO2*AvcaT+dSiO2)
              write(1012,*) runtime*tau/(year*1000),  AvIaT, AvImaT
              write(391,*) runtime*tau/(year*1000), (ta(DMFnd)*(Tl-Tsu)+Tsu)-273.15
              
              AmountS=int((AmountLeave+1)*MIDRDM)
              AmountI=AmountS !int((AmountS*(1.0-Uef)))
              AmountE=0 !int((AmountS*Uef))
              
              !if (AmountS.gt.(AmountI+AmountE)) then
               !   AmountI=AmountI+1
              !end if
              
              AILower=int(AmountI/MIDRDM)
              AELower=0!int(AmountE/MIDRDM)
              
              if ((AmountLeave+1).gt.(AILower+AELower)) then
                  AILower=AILower+1
              end if
              
              write(105,*) runtime*tau/(year*1000), (AmountI)*dz*delta/1000, ((AmountI)*dz*delta/1000)*((((MIDdleD/2)/1000)**2)*PI_8)/MIDRDM
              !write(106,*) runtime*tau/(year*1000), (AmountE)*dz*delta/1000, ((AmountE)*dz*delta/1000)*((((MIDdleD/2)/1000)**2)*PI_8)/MIDRDM
              
              write(107,*) runtime*tau/(year*1000), 'AS', AmountLeave+1
              write(107,*) runtime*tau/(year*1000), 'AI', AILower
              write(107,*) runtime*tau/(year*1000), 'AE', AELower
 
              write(110,*) runtime*tau/(year*1000), 'AS', AmountS
              write(110,*) runtime*tau/(year*1000), 'AI', AmountI
              write(110,*) runtime*tau/(year*1000), 'AE', AmountE
              
              write(100,*) 'Evacuated magma intruded at ', DMFnd*dz*delta/1000,'km'
              

              
              call eruptandshiftmass(j1b,j1b-AmountLeave,z1,z2,oldL3,oldL4,DMFnd,AmountI,AILower,AELower,AvhT, AvTT,&
                    &AvPaT, AvCaT, AvcsT, AvclT,ha,ta,pa,ca,cc,csa,cla,ica,icc,icsa,icla,ima,imc,imsa,imla,&
                    &hb,tb,pb,cb,csb,clb,icb,icsb,iclb,imb,imsb,imlb,clostim,density,&
                    & PartitionA,PartitionB,&
                    &runtime,year,tau,outcc,outcomp,baseflagP,topflagP, AvA1T, AvB1T,&
                    &AvC1T, AvTsT,A1,B1,C1,Ts,AvIaT,AvimaT, outci,outcim, outisot, outimt&
                    &,Tl,tsu,CloseTemp,upperD,midD,OGMIDLIM, OGUPLIM,d21,DMFNL,melt_seg_flag,melt_seg,pmag)
              
              
               ! pressureE=0
                pressureRTI=0
                !TotPressure=0
                taue=0
                DPL=5e20
                TRANSON=0
                t1=0
                t2=0
                hrtiP=0
                hlm1=0
                
                lambda=0
                lambdaP=0
                oldj3=5000000
                oldj4=5000000
                EruptTimeF=1

END IF     
    
    
    
    
! Find where the top and base of porosity is within the entire crust
    
102 do i=z1,z2
        if ((pb(i).gt.0).OR.(ca(i).ne.cc(i)).OR.(enth_nondim_to_dim(ha(i),Hs,Hl).ge.(cp*ts(i))).OR.((ta(i)*(Tl-Tsu)+Tsu).ge.Ts(i))) then
            if (i.gt.TopFlagP) then
                TopFlagP=i
            end if
            if (i.lt.BaseFlagP) then
                BaseFlagP=i
                
            end if
        end if
    end do
    
    !Update Velocity    
    if(velflag.gt.0)then
	    call velocitysolve(BaseFlagP,TopFlagP,dz,alpha,beta,pa,ws,wl,cla,max_melt_viscosity,min_melt_viscosity,cee,SLT,denscont,drho,aSio2,bSio2,cSio2,dSio2,melt_seg)
    else
    endif
        
      !! Start of eruption criteria
    
    !! check if magma can then erupt after a certain time frame 

if (EruptF.gt.0) then    
    if (runtime_erupt.ge.time_erupt_s) then
        top_Shall=5e6
        bottom_Shall=5e6
                    
        do i=z2,z1,-1
            if (i.gt.UpperDlim) then
                if ((pa(i).ge.cpa).AND.(pa(i+1).lt.cpa)) then
                    top_Shall=i
                    goto 651
                end if
            end if
        enddo
                    
    651 if (top_Shall.eq.5e6) then 
            runtime_erupt=0
            EruptTimeF=0
            write(100,*) 'There is no magma with a melt fraction greater than the CMF for eruption'
            goto 1022
        end if
                    
                    
        do i=top_Shall,z1,-1
            if ((pa(i).ge.cpa).AND.(pa(i-1).lt.cpa)) then
                    bottom_Shall=i
                    goto 751
            end if
        end do
                
                        
    751 if (bottom_Shall.eq.5e6) then
            runtime_erupt=0
            EruptTimeF=0
            write(100,*) 'There is no magma with a melt fraction greater than the CMF for eruption'
            goto 1022
        endif 
                    
                        
        write(100,*) 'Magma erupted at ', runtime*tau/(year*1000), 'ka' 
                    
        !!!!!!!!!!!!!!!!!!!!!!!! ERUPTION CRITIERA
                        !! calculate the averages of what erupts
        !1) find average comp and enth
                OcaE=0
                OccE=0
                OhE=0
                Oc1E=0
                
                do i=bottom_Shall,top_Shall
                    OcaE=OcaE+ca(i)
                    OhE=OhE+ha(i)
                    OccE=OccE+cc(i)
                    Oc1E=Oc1E+c1(i)
                end do 
                
                AvcaE=OcaE/(abs(top_Shall-bottom_Shall)+1)
                AvcdE=OccE/(abs(top_Shall-bottom_Shall)+1)
                AvhE=OhE/(abs(top_Shall-bottom_Shall)+1)
                Avc1E=Oc1E/(abs(top_Shall-bottom_Shall)+1)
                
                              
                if ((abs(avc1e-c1u)).lt.(abs(avc1e-c1l))) then
                        AvA1e=A1U
                        AvB1e=B1U
                        AvC1e=C1U
                        Avtse=TsU
                else
                      AvA1e=A1L
                      AvB1e=B1L
                      AvC1e=C1L
                      Avtse=TsL
                end if
                
            !2) find porosity
               call porosshift(AvhE,AvpaE,AvcaE,Hs,Hl,Lf,Cp,Tl,AvTse,Ae,AvA1e,AvB1e,AvC1e,A2,B2,C2)  
            
            !3) find composition of liquid and solid
               AvcsE=0
               AvclE=(1.0/AvpaE)*(AvcaE)
               if (AvclE.gt.1) then
                   AvclE=1.0
               end if
               
            !4) find average density of magma
               densMagLE=AvpaE*((1-ndsio(Avcl,aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(Avcl,aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-AvpaE)*((1-ndsio(Avcs,aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(Avcs,aSio2,bSio2,cSio2,dSio2)*DensSolB)      
               AvTE=AvA1e*AvclE**2+AvB1e*AvclE+AvC1e
               
               AvTE=temp_dim_to_nondim(AvTE,Tsu,Tl)   
               
        !! check it can erupt with density (a check)
               do i=top_Shall,z2
                   if (densMagLE.ge.density(i)) then
                       write(100,*) 'Density of the shallow magma is greater than the overlying crust - SHOULD WE ERUPT'
                   end if
               end do
               
                   
        !! write outputs of whats erupted (can use 106 for amounts - need a new one for comp)
       
               
               write(106,*) runtime*tau/(year*1000), abs(top_Shall-bottom_Shall)*dz*delta/1000, (abs(top_Shall-bottom_Shall)*dz*delta/1000)*((((MIDdleD/2)/1000)**2)*PI_8)/MIDRDM
               write(1060,*) runtime*tau/(year*1000),  aSio2-bSiO2*tanh(-cSiO2*AvcaE+dSiO2)
               
               if(outcountEM1.lt.10) then
		            write(filecountEM1,'(A1,A1,A1,A1,I1)') '0','0','0','0',outcountEM1
		            else if (outcountEM1.lt.100) then
		            write(filecountEM1,'(A1,A1,A1,I2)') '0','0','0',outcountEM1
		            else if (outcountEM1.lt.1000) then
		            write(filecountEM1,'(A1,A1,I3)') '0','0',outcountEM1
                    else if (outcountEM1.lt.10000) then
		            write(filecountEM1,'(A1,I4)') '0',outcountEM1
		            else
		            write(filecountEM1,'(I5)') outcountEM1
		            end if
		            outcountEM1 = outcountEM1 + 1
     
               
                    open(1061,file='Averages_Erupted'//filecountEM1//'.txt',status='old',err=21)
	                close(1061,status='delete')
21                  open(1061,file='Averages_Erupted'//filecountEM1//'.txt',status='new')      
                    write(1061,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
                    write(1061,*) 'Average density'
                    write(1061,*) DensMagLE
                    write(1061,*) 'Average comp'
                    write(1061,*) AvcaE
                    write(1061,*) 'Average temp'
                    write(1061,*) (AvTE*(Tl-Tsu)+Tsu)-273.15
                    write(1061,*) 'Average porosity'
                    write(1061,*) AvpaE
                    write(1061,*) 'Average cl'
                    write(1061,*) AvclE
                    close(1061)
                    
                    
                     !! erupt it!
                    
                     call eruptandshiftmass(top_Shall,bottom_Shall,z1,z2,oldL3,oldL4,DMFnd,0,0,abs(top_Shall-bottom_Shall),AvhE, AvTE,&
                    &AvPaE, AvCaE, AvcsE, AvclE,ha,ta,pa,ca,cc,csa,cla,ica,icc,icsa,icla,ima,imc,imsa,imla,&
                    &hb,tb,pb,cb,csb,clb,icb,icsb,iclb,imb,imsb,imlb,clostim,density,&
                    & PartitionA,PartitionB,&
                    &runtime,year,tau,outcc,outcomp,baseflagP,topflagP, AvA1T, AvB1T,&
                    &AvC1T, AvTsT,A1,B1,C1,Ts,AvIaT,AvimaT, outci,outcim, outisot, outimt&
                    &,Tl,tsu,CloseTemp,upperD,midD,OGMIDLIM, OGUPLIM,d21,DMFNL,melt_seg_flag,melt_seg,pmag)
                     
                     EruptTimeF=0
                     runtime_erupt=0
                    
                    
    endif
    
else
    EruptTimeF=0
    runtime_erupt=0
end if

    
    !! end of eruption criteria!!!!
    

     
     
     
      !Update old parameters and identify some key parameters.
     
1022     Transon=0   
    !Update old parameters and identify some key parameters.
    maxpi = 0
    maxupi = 0
    maxmpi=0
    maxlpi = 0    
    maxti = 0
    maxlp = 0
    maxmp=0
    maxup = 0
    maxlp = 0
    maxti = 0
    maxp = 0
    maxT = 0
    maxvel = 0
    maxdiff = 0
    maxdi = 0
    maxvi = 0
    
    midDlim=midD-xaboven
    if (midD.le.(z1+xaboven)) then
        midDlim=z1
    end if
    
    upperDlim=upperD-xaboven
    if (upperD.le.(z1+xaboven)) then
        upperDlim=z1
    end if
    

	do i = z1,z2
!! Find closure time which is when temp crosses the closure temp so T(t-1)<Tc and T(t)>Tc
        if(((ta(i)*(Tl-Tsu)+Tsu).ge.CloseTemp).and.((tb(i)*(Tl-Tsu)+Tsu).lt.CloseTemp))then
            clostim(i) = runtime*tau/year/1E6
        endif
        hb(i) = ha(i)
		tb(i) = ta(i)
		pb(i) = pa(i)
		cb(i) = ca(i)
		csb(i) = csa(i)
		clb(i) = cla(i)
		icb(i) = ica(i)
		iclb(i) = icla(i)
		icsb(i) = icsa(i)
		imb(i) = ima(i)
		imlb(i) = imla(i)
		imsb(i) = imsa(i)
		density(i) = pa(i)*((1-ndsio(cla(i),aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(cla(i),aSio2,bSio2,cSio2,dSio2)*DensLiqB)+(1-pa(i))*((1-ndsio(csa(i),aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(csa(i),aSio2,bSio2,cSio2,dSio2)*DensSolB)
		if((pa(i).gt.0.00001).AND.(pa(i).lt.1.0))then
		    denscont(i) = ((1-ndsio(csa(i),aSio2,bSio2,cSio2,dSio2))*DensSolA + ndsio(csa(i),aSio2,bSio2,cSio2,dSio2)*DensSolB) - ((1-ndsio(cla(i),aSio2,bSio2,cSio2,dSio2))*DensLiqA + ndsio(cla(i),aSio2,bSio2,cSio2,dSio2)*DensLiqB)
		else
		    denscont(i) = 0
		endif

! Find maximum porosity
        if(pa(i).gt.maxp)then
            maxp = pa(i)
            maxpi = i
         endif
! Find maximum temperature
        if(ta(i).lt.maxt)then
            maxt = ta(i)
            maxti = i
         endif
! Find maximum velocity
! NOTE THIS IS PORO-WEIGHTE
        if(abs((1-pb(i))*ws(i)).gt.maxvel)then
            maxvel = abs((1-pb(i))*ws(i))
            maxvi = i
        endif
         
    end do



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Output Data
!Only output time data after a certain time, and at a certain frequency
!if((runtime*tau/year/1000000).gt.writetime)then
    !if(timeoutcount.eq.timeoutfreq)then
!Temporal data at maximum porosity
   ! if(pa(maxpi).lt.0.0001)then
  !    write(50,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
 !   &  F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
!	&  runtime*tau/year/1000000,' ',dz*maxti*delta/1000,' ',(ta(maxti)*(Tl-Ts)+Ts)-273.5,' ',pa(maxti),' ',ca(maxti),' ',csa(maxti),' ',cla(maxti),&
!	&  ' ',aSio2-bSiO2*tanh(-cSiO2*ca(maxti)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(maxti)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(maxti)+dSiO2),&
	!&  ' ',ica(maxti),' ',icsa(maxti),' ',icla(maxti),' ',ima(maxti),' ',imsa(maxti),' ',imla(maxti),' ',ica(maxti)/ima(maxti),' ',icla(maxti)/imla(maxti),& 
	!&  ' ',clostim(maxti)    
   ! else
  !    write(50,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
 !   &  F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
!	&  runtime*tau/year/1000000,' ',dz*maxti*delta/1000,' ',(ta(maxpi)*(Tl-Ts)+Ts)-273.15,' ',pa(maxpi),' ',ca(maxpi),' ',csa(maxpi),' ',cla(maxpi),&
!	&  ' ',aSio2-bSiO2*tanh(-cSiO2*ca(maxpi)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(maxpi)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(maxpi)+dSiO2),&
!	&  ' ',ica(maxpi),' ',icsa(maxpi),' ',icla(maxpi),' ',ima(maxpi),' ',imsa(maxpi),' ',imla(maxpi),' ',ica(maxpi)/ima(maxpi),' ',icla(maxpi)/imla(maxpi),& 
!	&  ' ',clostim(maxpi) 
!    endif
!Temporal data at maximum porosity > SLT    
   ! if(pa(maxpi).gt.SLT)then
  !  write(60,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
 !   &  F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
!	&  runtime*tau/year/1000000,' ',dz*maxti*delta/1000,' ',(ta(maxpi)*(Tl-Ts)+Ts)-273.15,' ',pa(maxpi),' ',ca(maxpi),' ',csa(maxpi),' ',cla(maxpi),&
!	&  ' ',aSio2-bSiO2*tanh(-cSiO2*ca(maxpi)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(maxpi)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(maxpi)+dSiO2),&
!	&  ' ',ica(maxpi),' ',icsa(maxpi),' ',icla(maxpi),' ',ima(maxpi),' ',imsa(maxpi),' ',imla(maxpi),' ',ica(maxpi)/ima(maxpi),' ',icla(maxpi)/imla(maxpi),&
!	&  ' ',clostim(maxpi)  
!    endif
!Temporal data at all porosity > SLT     
!	do i = z1,z2
   ! if(pa(i).gt.SLT)then
  !  write(70,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
 !   &  F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
!	&  runtime*tau/year/1000000,' ',dz*i*delta/1000,' ',(ta(i)*(Tl-Ts)+Ts)-273.15,' ',pa(i),' ',ca(i),' ',csa(i),' ',cla(i),&
!	&  ' ',aSio2-bSiO2*tanh(-cSiO2*ca(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(i)+dSiO2),&
!	&  ' ',ica(i),' ',icsa(i),' ',icla(i),' ',ima(i),' ',imsa(i),' ',imla(i),' ',ica(i)/ima(i),' ',icla(i)/imla(i),' ',clostim(i)  
!	freqcount = freqcount + 1
!    endif
!	end do
!Temporal data at specific depth
!write(*,*)depthi,dz*depthi*delta/1000
  !  i=depthi
   ! if(depthoutput.gt.0)then
  !  write(80,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
 !   &  F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
!	&  runtime*tau/year/1000000,' ',(ta(i)*(Tl-Ts)+Ts)-273.15,' ',pa(i),' ',ca(i),' ',csa(i),' ',cla(i),&
!	&  ' ',aSio2-bSiO2*tanh(-cSiO2*ca(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(i)+dSiO2),&
!	&  ' ',ica(i),' ',icsa(i),' ',icla(i),' ',ima(i),' ',imsa(i),' ',imla(i),' ',ica(i)/ima(i),' ',icla(i)/imla(i),' ',clostim(i)  
!	endif
 !   timeoutcount = 0
!	else
!	timeoutcount=timeoutcount+1
!	endif
!else
!endif
if (maxpor.gt.0) then    
    do i=z1,midDlim_LOWERONLY
        if (pa(i).gt.maxlp) then
            maxlp=pa(i)
            maxlpi=i
        end if
    end do
    
    do i=midDlim, upperDlim-1
        if (pa(i).gt.maxmp) then
            maxmp=pa(i)
            maxmpi=i
        end if
    end do
    

    do i=upperDlim,z2
        if (pa(i).gt.maxup) then
            maxup=pa(i)
            maxupi=i
        end if
    end do



write(98,*) maxlpi, maxlp
write(981,*)maxmpi, maxmp
write(99,*) maxupi, maxup
end if

!Spatial data - data only written every n years specified at input
if (runtime/real(outcount).gt.OutRate-dtl/real(outcount)) then
	if (runtime/real(outcount).lt.OutRate+dtl/real(outcount)) then
		write(*,'(A12,I4)'), 'Data Output: ',outcount
		write(100,'(A12,I4)'), 'Data Output: ',outcount
		write(*,'(A1,F7.3,A10)'),' ',100.0*runtime/time,'% Complete'
		write(100,'(A1,F7.3,A10)'),' ',100.0*runtime/time,'% Complete'
		write(*,'(A50,F7.3,A2)')'Bulk mass balance error over interval of interest ',cerror,' %'
		write(100,'(A50,F17.13,A2)')'Bulk mass balance error over interval of interest ',cerror,' %'
	    write(*,'(A53,F7.3,A2)')'Isotope (crust) mass balance error over interval of interest ',ierror,' %'
        write(*,'(A53,F7.3,A2)')'Isotope (mantle) mass balance error over interval of interest ',imerror,' %'
		write(100,'(A53,F17.13,A2)')'Isotope mass balance error over interval of interest ',ierror,' %'
		write(*,'(A17,F17.13)')'Maximum velocity ',maxvel
		write(100,'(A17,F12.9)')'Maximum velocity ',maxvel
		write(*,'(A16,F17.13)')'Upwind weighting ',weighttemp
		write(100,'(A16,F17.13)')'Upwind weighting ',weighttemp
		write(*,'(A10,F17.13)')'Timestep ',dtl
		write(100,'(A10,F12.9)')'Timestep ',dtl
        
        write(103,*) runtime*tau/(year*1000),cerror
        write(1031,*) runtime*tau/(year*1000),ierror
        write(1032,*) runtime*tau/(year*1000),imerror
        
        write(108,*) runtime*tau/(year*1000), TopFlagP, BaseFlagP
        
        write(119,*) runtime*tau/(year*1000), pressureE
        write(120,*) runtime*tau/(year*1000), pressureRTI
        write(121,*) runtime*tau/(year*1000), Totpressure
        write(122,*) runtime*tau/(year*1000), DPL_OG
        write(123,*) runtime*tau/(year*1000), pressureEL
        write(124,*) runtime*tau/(year*1000), pressureRTIL
        write(125,*) runtime*tau/(year*1000), TotpressureL
        write(126,*) runtime*tau/(year*1000), LDPL_OG
        write(129,*) runtime*tau/(year*1000), t1L/(year*1000)
        write(130,*) runtime*tau/(year*1000), t1/(year*1000)
        write(135,*) runtime*tau/(year*1000), denRL-LdensMagL
        write(136,*) runtime*tau/(year*1000), denR-densMagL
        write(137,*) runtime*tau/(year*1000), hRTIL
        write(138,*) runtime*tau/(year*1000), hRTI
        write(139,*) runtime*tau/(year*1000), tauEL
        write(140,*) runtime*tau/(year*1000), tauE
        write(141,*) runtime*tau/(year*1000), t2L/(year*1000)
        write(142,*) runtime*tau/(year*1000), t2/(year*1000)
        write(143,*) runtime*tau/(year*1000), MagViscL, AvporBL
        write(144,*) runtime*tau/(year*1000), MagVisc, AvporB
        
        write(319,*) runtime*tau/(year*1000), upperDlim
        write(320,*) runtime*tau/(year*1000), midDlim
        write(3201,*) runtime*tau/(year*1000), midDlim_LOWERONLY
        
        write(321,*) runtime*tau/(year*1000),oldj3, oldj4
        write(322,*)runtime*tau/(year*1000),j3, j4
        write(323,*)runtime*tau/(year*1000),oldL3, oldL4
        write(324,*)runtime*tau/(year*1000),l3,l4
        
        write(325,*)runtime*tau/(year*1000),lowd
        write(326,*)runtime*tau/(year*1000),middled
        write(327,*)runtime*tau/(year*1000),lowcrustm
        write(328,*)runtime*tau/(year*1000),midcrustm
        write(329,*)runtime*tau/(year*1000),oguplim
        write(330,*)runtime*tau/(year*1000),ogmidlim
        
		if(outcount.lt.10) then
		write(filecount,'(A1,A1,A1,I1)') '0','0','0',outcount
		else if (outcount.lt.100) then
		write(filecount,'(A1,A1,I2)') '0','0',outcount
		else if (outcount.lt.1000) then
		write(filecount,'(A1,I3)') '0',outcount
		else
		write(filecount,'(I4)') outcount
		end if
		outcount = outcount + 1
	else
		goto 500
	end if
else
	goto 500
end if
	open(1000,file='output'//filecount//'.txt',status='old',err=10)
	close(1000,status='delete')
10	open(1000,file='output'//filecount//'.txt',status='new')

	write(1000,'(A4,F10.5,A3)'), 'Time', runtime*tau/(year*1000),' ka'
	write(1000,'(26A20)') 'Distance (km)','Enthalpy (J)','Temperature (C)','Porosity','Melt Seg','P mag','Solid Velocity (m/s)','Liquid Velocity (m/s)','Bulk Composition','Solid Composition','Liquid Composition',&
	&'Bulk SiO2 (%)','Solid SiO2 (%)','Liquid SiO2 (%)',&
	&'Bulk Crust I','Solid Crust I','Liquid Crust I','Bulk Mantle I','Solid Mantle I','Liquid Mantle I', 'Bulk IC/BulkIM', 'Liquid IC/LiquidIM',&
    &'ts','Closure Time (My)','Density (kg/m3)','Den contrast (kg/m3)'
    do i = z1,z2
!        if((ta(i)*(Tl-Ts)+Ts).gt.CloseTemp)then
!	write(1000,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
!	&F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A6)') &
!	&dz*i*delta/1000,' ',ha(i),' ',(ta(i)*(Tl-Ts)+Ts)-273.15,' ',pa(i),' ',ws(i),' ',wl(i),' ',ca(i),' ',csa(i),' ',cla(i),&
!	&' ',aSio2-bSiO2*tanh(-cSiO2*ca(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(i)+dSiO2),&
!	&' ',ica(i),' ',icsa(i),' ',icla(i),' ',ima(i),' ',imsa(i),' ',imla(i),' Undef'
!       else
	write(1000,'(F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,&
	&F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10,A1,F19.10)') &
	&dz*i*delta/1000,' ',ha(i),' ',(ta(i)*(Tl-Tsu)+Tsu)-273.15,' ',pa(i),' ',melt_seg(i),' ',pmag(i),' ',ws(i),' ',wl(i),' ',ca(i),' ',csa(i),' ',cla(i),&
    &' ',aSio2-bSiO2*tanh(-cSiO2*ca(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*csa(i)+dSiO2),' ',aSio2-bSiO2*tanh(-cSiO2*cla(i)+dSiO2),&
	&' ',ica(i),' ',icsa(i),' ',icla(i),' ',ima(i),' ',imsa(i),' ',imla(i),' ',ica(i)/ima(i),' ',icla(i)/imla(i),' ',ts(i),' ',clostim(i),' ',density(i),' ',denscont(i)
!	    endif
	end do
    close(1000)
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

500 dtadv = dz/(tmultad1*maxvel)
    dtartdiff = (dz**2)/(tmultad2*maxdiff*maxvel)
    dtl = min(dtadv,dtartdiff,dtdiff)
    
    !dtl = (dz**2)/((tmultdiff*keff)+(tmultad1*maxvel*dz)+(tmultad2*maxvel*maxdiff))
!    dtadv = dz/(tmult*maxvel)
!    dtartdiff = (dz**2)/(tmult*maxdiff*maxvel)
!    dtl = min(dtadv,dtartdiff,dtdiff)
!    write(100,*)dtl,dtadv,dtartdiff,dtdiff,maxdiff
!    if(dtadv.lt.dtdiff)then
!        dtl = dtadv
!    else
!        dtl = dtdiff
!    endif
    runtime = runtime + dtl
    runtimeT = runtimeT + dtl
    dsilltime = dsilltime + dtl
    
    if (EruptTimeF.gt.0) then
        runtime_erupt=runtime_erupt+dtl
    end if
    
    
    
	j = j+1

end do



write(*,*)'Total number of records in frequency data file: ',freqcount
write(100,*)'Total number of records in frequency data file: ',freqcount

close(50)
close(60)
close(70)
close(80)
close(100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end program

! SUBROUTINES AND FUNCTIONS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sillemplace(z1,z2,d1,d2,SillNodes,Accrete,&
					  &ha,ta,pa,ca,csa,cla,ica,icsa,icla,ima,imsa,imla,&
					  &hb,tb,pb,ws,wl,cb,cc,csb,clb,icb,icc,imc,icsb,iclb,imb,&
                      &imsb,imlb,SillEnth,SillTemp,SillPhi,SillComp,SillIso1,&
					  &SillIso2,Ts,Tl,PartitionA,PartitionB,Ae,A1,B1,C1,A2,B2&
                    &,C2,dz,delta,clostim,SN,A1U,B1U,C1U,Tsu,runtime,tau,year,CloseTemp, melt_seg,pmag)
implicit none
integer i,n,z1,z2,d1,d2,sillnodes,silln,sillcount,accrete,SN
parameter (n=1000000)

real(8) hb(-n:n)	    !Enthalpy
real(8) tb(-n:n)	    !Temperature
real(8) pb(-n:n)	    !Porosity
real(8) ws(-n:n),wl(-n:n)	!Velocity
real(8) cb(-n:n)	    !Bulk Composition
real(8) cc(-n:n)	    !Bulk Composition Check
real(8) csb(-n:n)	    !Solid Composition
real(8) clb(-n:n)	    !Liquid Composition
real(8) icb(-n:n)	    !Isotope Bulk Composition
real(8) icc(-n:n)	    !Isotope Bulk Composition Check
real(8) icsb(-n:n)	    !Solid Composition
real(8) iclb(-n:n)	    !Liquid Composition
real(8) imb(-n:n)	    !Isotope Bulk Composition
real(8) imsb(-n:n)	    !Solid Composition
real(8) imlb(-n:n)	    !Liquid Composition
real(8) imc(-n:n)       ! isotope check
real(8) ha(-n:n)	    !Enthalpy
real(8) ta(-n:n)	    !Temperature
real(8) pa(-n:n)	    !Porosity
real(8) ca(-n:n)	    !Bulk Composition
real(8) csa(-n:n)	    !Solid Composition
real(8) cla(-n:n)	    !Liquid Composition
real(8) ica(-n:n)	    !Isotope Bulk Composition
real(8) icsa(-n:n)	    !Solid Composition
real(8) icla(-n:n)	    !Liquid Composition
real(8) ima(-n:n)	    !Isotope Bulk Composition
real(8) imsa(-n:n)	    !Solid Composition
real(8) imla(-n:n)	    !Liquid Composition
real(8) clostim(-n:n)	!Closure time
real(8) A1(-n:n)
real(8) B1(-n:n)
real(8) C1(-n:n)
real(8) Ts(-n:n)
real(8) melt_seg(-n:n)
real(8) melt_seg1(-n:n)
real(8) pmag(-n:n)

real(8) SillEnth,SillTemp,SillPhi,SillComp,SillIso1,SillIso2,dz,delta

real(8) Tl,Ae,A2,B2,C2
real(8) A1U, B1U, C1U, TsU
real(8) PartitionA,PartitionB,shiftn

real(8) SolidComp,LiquidComp
real(8) IsoSolid,IsoLiquid
real(8) temp_nondim_to_dim

real(8) runtime, year, tau, closeTemp

!write(*,*)z1,z2,d1,d2
!Over-Accretion
!if (accrete.eq.0) then

shiftn=sillnodes

	write(*,*) 'Sill emplaced at: ', d2*dz*delta/1000,' km'
	write(100,*) 'Sill emplaced at: ', d2*dz*delta/1000,' km'
    write(*,*) d2*dz*delta/1000
	do i = z1,d2
		hb(i-shiftn) = hb(i)
		tb(i-shiftn) = tb(i)
		pb(i-shiftn) = pb(i)
		
		ws(i-shiftn) = ws(i)
		wl(i-shiftn) = wl(i)
		
		cb(i-shiftn) = cb(i)
        cc(i-shiftn) = cc(i)
		csb(i-shiftn) = csb(i)
		clb(i-shiftn) = clb(i)
		
		icb(i-shiftn) = icb(i)
		icc(i-shiftn) = icc(i)
		icsb(i-shiftn) = icsb(i)
		iclb(i-shiftn) = iclb(i)
		
		imb(i-shiftn) = imb(i)
		imsb(i-shiftn) = imsb(i)
		imlb(i-shiftn) = imlb(i)
		imc(i-shiftn) = imc(i)
        
		ha(i-shiftn) = ha(i)
		ta(i-shiftn) = ta(i)
		pa(i-shiftn) = pa(i)
		
		ca(i-shiftn) = ca(i)
		csa(i-shiftn) = csa(i)
		cla(i-shiftn) = cla(i)
		
		ica(i-shiftn) = ica(i)
		icsa(i-shiftn) = icsa(i)
		icla(i-shiftn) = icla(i)
        

		ima(i-shiftn) = ima(i)
		imsa(i-shiftn) = imsa(i)
		imla(i-shiftn) = imla(i)
        
        pmag(i-shiftn) = pmag(i)
		
        clostim(i-shiftn) = clostim(i)
        
        A1(i-shiftn)=A1(i)
        B1(i-shiftn)=B1(i)
        C1(i-shiftn)=C1(i)
        Ts(i-shiftn)=Ts(i)
        
        Melt_seg(i-shiftn)=Melt_seg(i)

	end do
	do i = d2-(sillnodes-1),d2
		hb(i) = SillEnth
		tb(i) = SillTemp
		pb(i) = SillPhi

		ws(i) = 0
		wl(i) = 0

		cb(i) = SillComp
		cc(i) = SillComp
		csb(i) = SolidComp(temp_nondim_to_dim(SillTemp,Tsu,Tl),SillComp,Ae,A1(i),B1(i),C1(i),A2,B2,C2,Ts(i))
		clb(i) = LiquidComp(temp_nondim_to_dim(SillTemp,Tsu,Tl),SillComp,Ts(i),Ae,A1(i),B1(i),C1(i),A2,B2,C2)

		icb(i) = SillIso1
		icc(i) = SillIso1
		icsb(i) = IsoSolid(pb(i),icb(i),cb(i),clb(i),csb(i),PartitionA,PartitionB)		
		iclb(i) = IsoLiquid(pb(i),icb(i),cb(i),clb(i),csb(i),PartitionA,PartitionB)
		imb(i) = SillIso2
        imc(i) = SillIso2
		imsb(i) = IsoSolid(pb(i),imb(i),cb(i),clb(i),csa(i),PartitionA,PartitionB)	
		imlb(i) = IsoLiquid(pb(i),imb(i),cb(i),clb(i),csa(i),PartitionA,PartitionB)
        pmag(i) = 1
		
		if(((SillTemp*(Tl-Tsu)+Tsu).ge.CloseTemp))then
            clostim(i) = runtime*tau/year/1E6
        else
            clostim(i)=-999
        end if
        
        A1(i)=A1U
        B1(i)=B1U
        C1(i)=C1U
        Ts(i)=TsU
        Melt_seg(i)=1
    end do
    write(*,*) z1-sillnodes, d2-sillnodes, d2-(sillnodes-1),d2
	z1 = z1-shiftn
    SN = sillnodes
!	d1 = d1-sillnodes ! THIS WAS COMMENTED OUT INCLUDE IF FAILS

!Under-Accretion
!else if (accrete.eq.1) then
!	write(*,*) 'Sill emplaced at: ',d2-sillcount*sillnodes	
!	do i = z1,(d2-(sillcount-1)*sillnodes)
	!	ha(i-sillnode) = ha(i)
	!	hb(i-sillnode) = hb(i)
	!	ta(i-sillnode) = ta(i)
	!	tb(i-sillnode) = tb(i)
	!	pa(i-sillnode) = pa(i)
	!	pb(i-sillnode) = pb(i)
	!	ws(i-sillnode) = ws(i)
	!	wl(i-sillnode) = wl(i)
	!end do
	!do i = (2-sillcount*sillnodes + 1),(d2-(sillcount-1)*sillnodes)
	!	ha(i) = SillEnth
	!	hb(i) = SillEnth
	!	ta(i) = SillTemp
	!	tb(i) = SillTemp
	!	pa(i) = SillPhi
	!	pb(i) = SillPhi
	!	ws(i) = 0
	!	wl(i) = 0
	!end do
	!z1 = z1-sillnodes
!end if

end subroutine
!========================================================
subroutine enthsolve(z1,z2,dt,dz,keff,stefan,ha,hb,tb,pb,ws,wl)
implicit none

integer i,n,z1,z2
parameter (n=1000000)
real(8) dt,dz
real(8) keff,stefan
real(8) ha(-n:n),hb(-n:n),tb(-n:n),pa(-n:n),pb(-n:n),ws(-n:n),wl(-n:n)
real(8) conduct(-n:n),advect(-n:n)

! BC on enthalpy/temp need to be considered carefully

do i = z1+1,z2-1
!	if (i.eq.z2) then
! Zero gradient in enthalpy at boundary
!		conduct(i) = keff*(dt/dz**2.0)*(-tb(i) + tb(i-1))
		!advect(i) =  stefan*(dt/dz)*(-(1.0-pb(i))*ws(i))
!		advect(i) = 0
!	elseif (i.eq.z1) then
!		conduct(i) = keff*(dt/dz**2.0)*(tb(i+1) - 2.0*tb(i) + tb(i-1))
!		conduct(i) = keff*(dt/dz**2.0)*(tb(i+1) - tb(i))

!		advect(i) = stefan*(dt/dz)*((1.0-pb(i+1))*ws(i+1)-(1.0-pb(i))*ws(i))
!	else
	conduct(i) = keff*(dt/dz**2.0)*(tb(i+1) - 2.0*tb(i) + tb(i-1))	
	advect(i) = stefan*(dt/(2*dz))*((1.0-pb(i+1))*ws(i+1)-(1.0-pb(i-1))*ws(i-1))
!	end if

	ha(i) = hb(i) + conduct(i) + advect(i)

end do

! Fixed enthalpy/temp at top and base - look at McKenzie?  
	ha(z1) = hb(z1)
	ha(z2) = hb(z2)

end subroutine
!========================================================
!========================================================
!COMPSOLVE
subroutine compsolve(BaseFlagP, TopFlagP, z1,z2,dt,dz,pa,pb,ca,cb,ws,wl,cs,cl,maxvel,maxdiff,cc,cerror,cerrorU,cerrorL,baseflagU,topflagU,baseflagL,topflagL,outcomp,outcc,&
&diffexp,diffmax,diffmin,ia,ib,il,is,ic,ierror,outisot,outci,ima,imb,iml,ims,imc,imerror,outimt,outcim,CompChoice)
implicit none

integer i,n
parameter (n=1000000)
integer z1,z2,baseflagU,topflagU,baseflagL,topflagL,T , baseflagP, topflagP

real(8) dt,dz,maxvel,maxdiff,cerror,centsum,ccsum,upsum,downsum,centerror,uperror,downerror,cerrorU,cerrorL
!real(8) cmin, cmaxm, laxdum, laxdumexth, laxdumextl

real(8) pa(-n:n),pb(-n:n),ws(-n:n),wl(-n:n)
real(8) ca(-n:n),cb(-n:n)
real(8) cl(-n:n),cs(-n:n),cc(-n:n)
real(8) upc(-n:n),downc(-n:n),centc(-n:n)
real(8) ia1(-n:n), ima1(-n:n)
real(8) outcomp, outcc

!! Define isotopes
real(8) ia(-n:n),ib(-n:n)
real(8) ima(-n:n), imb(-n:n)
real(8) il(-n:n),is(-n:n),ic(-n:n)
real(8) iml(-n:n), ims(-n:n), imc(-n:n)
real(8) iup(-n:n), imup(-n:n)
real(8) ierror,iccsum,iupsum,iuperror
real(8) imerror, imcsum, imupsum, imuperror
real(8) outisot,outci
real(8) outimt, outcim
real(8) ia1sum, iaerror, ia1error
real(8) ima1sum, imaerror, ima1error

!! Calculate the diffusion terms locally instead of in main body and define variables here
real(8) cdiffdum, idiffdum,imdiffdum,diffexp,diffmax,diffmin

real(8) CompChoice

!write(*,*)laxm, laxe
T=0


if (CompChoice.eq.0) then 
!!! old solver - composition is found first, then isotopes using the method with the 
    !best mass conservation. 
    
    do i=BaseFlagP, TopFlagP
        if (pb(i).gt.0) then 
            
                cdiffdum = (((2.0*cb(i)-1.0)**diffexp)*(diffmax-diffmin))+diffmin
        
                centc(i) = cb(i)&
	                        &-(dt/(2*dz))*(cs(i+1)*ws(i+1)*(1.0-pb(i+1)) - cs(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                    &+(dt/(2*dz))*(cl(i+1)*ws(i+1)*(1.0-pb(i+1)) - cl(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                    !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                            &+cdiffdum*maxvel*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
        !        
                upc(i) = cb(i)&
	                        &-(dt/(dz))*(cs(i+1)*ws(i+1)*(1.0-pb(i+1)) - cs(i)*ws(i)*(1.0-pb(i)))&
		                    &+(dt/(dz))*(cl(i+1)*ws(i+1)*(1.0-pb(i+1)) - cl(i)*ws(i)*(1.0-pb(i)))&
                    !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                            &+cdiffdum*maxvel*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
        !        
   	            downc(i) = cb(i)&
	                        &-(dt/(dz))*(cs(i)*ws(i)*(1.0-pb(i)) - cs(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                    &+(dt/(dz))*(cl(i)*ws(i)*(1.0-pb(i)) - cl(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                    !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                            &+cdiffdum*maxvel*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
        else
            centc(i)=cb(i)
            downc(i)=cb(i)
            upc(i)=cb(i)
        end if
        

    ! This term violates conservation of mass but ensures a few nodes that are out of bounds do not destroy solution
        if(centc(i).lt.0.00001) then
	        !write(*,*)'Warning Cb (centered) below zero at node ',i,' with value ',centc(i)
	        !write(100,*)'Warning Cb (centered) below zero at node ',i,' with value ',centc(i)
	        centc(i) = 0.00001
        else if(centc(i).gt.0.99999) then
            !write(*,*)'Warning Cb (centered) above one at node ',i,' with value ',centc(i)
            !write(100,*)'Warning Cb (centered) above one at node ',i,' with value ',centc(i)
	        centc(i) = 0.99999
        end if
	
        if(upc(i).lt.0.00001) then
	        !write(*,*)'Warning Cb (upwind) below zero at node ',i,' with value ',upc(i)
        !		write(100,*)'Warning Cb (upwind) below zero at node ',i,' with value ',upc(i)
	        upc(i) = 0.00001
        else if(upc(i).gt.0.99999) then
            !write(*,*)'Warning Cb (upwind) above one at node ',i,' with value ',upc(i)
            !write(100,*)'Warning Cb (upwind) above one at node ',i,' with value ',upc(i)
	        upc(i) = 0.99999
        end if
	
        if(downc(i).lt.0.00001) then
	        !write(*,*)'Warning Cb (downwind) below zero at node ',i,' with value ',downc(i)
	        !write(100,*)'Warning Cb (downwind) below zero at node ',i,' with value ',downc(i)
	        downc(i) = 0.00001
        else if(downc(i).gt.0.99999) then
            !write(*,*)'Warning Cb (downwind) above one at node ',i,' with value ',downc(i)
            !write(100,*)'Warning Cb (downwind) above one at node ',i,' with value ',downc(i)
	        downc(i) = 0.99999
        end if        
    end do
! Check mass balance - overall
    centsum=0
    upsum=0
    downsum=0
    ccsum=0
    centerror=0
    uperror=0
    downerror=0    
    

    do i = BaseFlagP,TopFlagP
        centsum = centsum + centc(i)
        upsum = upsum + upc(i) 
        downsum = downsum + downc(i) 
        ccsum = ccsum + cc(i)
    end do
    
    centerror = (centsum+outcomp-(ccsum+outcc))/(ccsum+outcc)
    uperror = (upsum+outcomp-(ccsum+outcc))/(ccsum+outcc)
    downerror = (downsum+outcomp-(ccsum+outcc))/(ccsum+outcc)  
    
    if(abs(centerror).gt.abs(uperror))then
    if(abs(uperror).gt.abs(downerror))then
        do i = BaseFlagP,TopFlagP
            if (pb(i).gt.0) then
                ca(i)=downc(i)/(1+downerror)
            else
                ca(i)=downc(i)
            end if
            
        enddo
        cerror=downerror*100      
        T=1
    else
        do i = BaseFlagP,TopFlagP
            if (pb(i).gt.0) then
                ca(i)=upc(i)/(1+uperror)
            else
                ca(i)=upc(i)
            end if
            
        enddo
        cerror=uperror*100 
        T=2
    endif
    else if(abs(centerror).gt.abs(downerror))then 
        do i = BaseFlagP,TopFlagP
            if (pb(i).gt.0) then 
                ca(i)=downc(i)/(1+downerror)
            else
                ca(i)=downc(i)
            end if
            
        enddo
        cerror=downerror*100
        T=1
    else
        do i = BaseFlagP,TopFlagP
            if (pb(i).gt.0) then
                ca(i)=centc(i)/(1.0+centerror)
            else
                ca(i)=centc(i)
            end if
        enddo
        cerror=centerror*100
        T=3
    end if
    
   do i=BaseFlagP,TopFlagP
       if (pb(i).gt.0) then
           idiffdum = (diffmax*(diffmin**(diffexp/10.0))/((ib(i)+1.0-diffmin)**(diffexp/10.0)))+diffmin
           imdiffdum = (diffmax*(diffmin**(diffexp/10.0))/((imb(i)+1.0-diffmin)**(diffexp/10.0)))+diffmin
          
            if (T.eq.1) then
                !isotopes - down error
                !crust
                ia1(i)= ib(i)&
	                    &-(dt/(dz))*(is(i)*ws(i)*(1.0-pb(i)) - is(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                &+(dt/(dz))*(il(i)*ws(i)*(1.0-pb(i)) - il(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                        &+idiffdum*maxvel*(dt/dz**2.0)*(ib(i+1)-2.0*ib(i)+ib(i-1))
                
                !mantle
                ima1(i)= imb(i)&
	                    &-(dt/(dz))*(ims(i)*ws(i)*(1.0-pb(i)) - ims(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                &+(dt/(dz))*(iml(i)*ws(i)*(1.0-pb(i)) - iml(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                        &+imdiffdum*maxvel*(dt/dz**2.0)*(imb(i+1)-2.0*imb(i)+imb(i-1))
        
            else if (T.eq.2) then
                !isotopes - up error
                !crust
                ia1(i) = ib(i)&
	                    &-(dt/(dz))*(is(i+1)*ws(i+1)*(1.0-pb(i+1)) - is(i)*ws(i)*(1.0-pb(i)))&
		                &+(dt/(dz))*(il(i+1)*ws(i+1)*(1.0-pb(i+1)) - il(i)*ws(i)*(1.0-pb(i)))&
                !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                        &+idiffdum*maxvel*(dt/dz**2.0)*(ib(i+1)-2.0*ib(i)+ib(i-1))
                
                !mantle
                ima1(i) = imb(i)&
	                    &-(dt/(dz))*(ims(i+1)*ws(i+1)*(1.0-pb(i+1)) - ims(i)*ws(i)*(1.0-pb(i)))&
		                &+(dt/(dz))*(iml(i+1)*ws(i+1)*(1.0-pb(i+1)) - iml(i)*ws(i)*(1.0-pb(i)))&
                !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                        &+imdiffdum*maxvel*(dt/dz**2.0)*(imb(i+1)-2.0*imb(i)+imb(i-1))
    
            else if (T.eq.3) then
                !isotopes - cent error 
                !crust
                ia1(i) = ib(i)&
	                    &-(dt/(2*dz))*(is(i+1)*ws(i+1)*(1.0-pb(i+1)) - is(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                &+(dt/(2*dz))*(il(i+1)*ws(i+1)*(1.0-pb(i+1)) - il(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                        &+idiffdum*maxvel*(dt/dz**2.0)*(ib(i+1)-2.0*ib(i)+ib(i-1))
                !mantle
                ima1(i) = imb(i)&
	                    &-(dt/(2*dz))*(ims(i+1)*ws(i+1)*(1.0-pb(i+1)) - ims(i-1)*ws(i-1)*(1.0-pb(i-1)))&
		                &+(dt/(2*dz))*(iml(i+1)*ws(i+1)*(1.0-pb(i+1)) - iml(i-1)*ws(i-1)*(1.0-pb(i-1)))&
                !	    &+(dz*laxdum)*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
                        &+imdiffdum*maxvel*(dt/dz**2.0)*(imb(i+1)-2.0*imb(i)+imb(i-1))
            end if
       else
           ia1(i)=ib(i)
           ima1(i)=imb(i)
       end if
       
            
        if (ia1(i).lt.0.00001) then
	       ! write(*,*)'Warning Ib below zero at node ',i,' with value ',ia(i)
		    !write(100,*)'Warning Ib below zero at node ',i,' with value ',ia(i)
		    ia1(i) = 0.00001
        end if	
       
        if (ima1(i).lt.0.00001) then
	      !  write(*,*)'Warning Imb below zero at node ',i,' with value ',ima(i)
		    !write(100,*)'Warning Imb below zero at node ',i,' with value ',ima(i)
		    ima1(i) = 0.00001
	    end if
        
           
   end do
   
   ia1sum=0
   ia1error=0
   iccsum=0
   
   ima1sum=0
   ima1error=0
   imcsum=0
   
    do i = BaseFlagP,TopFlagP
        ia1sum = ia1sum + ia1(i)
        iccsum = iccsum + ic(i)
        
        ima1sum = ima1sum + ima1(i)
        imcsum = imcsum + imc(i)
    end do
  !  write(*,*) ia1sum, ima1sum
    
    
    ia1error = (ia1sum+outisot-(iccsum+outci))/(iccsum+outci)
    ima1error = (ima1sum+outimt-(imcsum+outcim))/(imcsum+outcim)
   
   !write(*,*) ima1sum, outimt, imcsum, outcim
    do i=BaseFlagP, TopFlagP
        if (pb(i).gt.0) then
            ia(i)=ia1(i)/(1.0+ia1error)
            ima(i)=ima1(i)/(1.0+ima1error)
        else
            ia(i)=ia1(i)
            ima(i)=ima1(i)
        end if
    end do

    ierror=ia1error*100
    imerror=ima1error*100
    

    else
        !! Matts new code
      
       do i = BaseFlagP, TopFlagP
           if (pb(i).gt.0) then
               
                ! Calculate artificial diffusion here
                 cdiffdum = (((2*cb(i)-1)**diffexp)*(diffmax-diffmin))+diffmin
                 idiffdum = (diffmax*(diffmin**(diffexp/10))/((ib(i)+1-diffmin)**(diffexp/10)))+diffmin  
                 imdiffdum = (diffmax*(diffmin**(diffexp/10))/((imb(i)+1-diffmin)**(diffexp/10)))+diffmin  

               ! Upwinding   
                upc(i) = cb(i)&
	            &-(dt/(dz))*(cs(i+1)*ws(i+1)*(1.0-pb(i+1)) - cs(i)*ws(i)*(1.0-pb(i)))&
		        &+(dt/(dz))*(cl(i+1)*ws(i+1)*(1.0-pb(i+1)) - cl(i)*ws(i)*(1.0-pb(i)))&
                &+cdiffdum*maxvel*(dt/dz**2.0)*(cb(i+1)-2.0*cb(i)+cb(i-1))
        
                iup(i) = ib(i)&
	            &-(dt/(dz))*(is(i+1)*ws(i+1)*(1.0-pb(i+1)) - is(i)*ws(i)*(1.0-pb(i)))&
		        &+(dt/(dz))*(il(i+1)*ws(i+1)*(1.0-pb(i+1)) - il(i)*ws(i)*(1.0-pb(i)))&
                &+idiffdum*maxvel*(dt/dz**2.0)*(ib(i+1)-2.0*ib(i)+ib(i-1))
                
                imup(i) = imb(i)&
	            &-(dt/(dz))*(ims(i+1)*ws(i+1)*(1.0-pb(i+1)) - ims(i)*ws(i)*(1.0-pb(i)))&
		        &+(dt/(dz))*(iml(i+1)*ws(i+1)*(1.0-pb(i+1)) - iml(i)*ws(i)*(1.0-pb(i)))&
                &+imdiffdum*maxvel*(dt/dz**2.0)*(imb(i+1)-2.0*imb(i)+imb(i-1))
       
           else 
               
                upc(i) = cb(i) 
                iup(i) = ib(i)
                imup(i) = imb(i)

            endif  
        
            if(upc(i).lt.0.00001) then
		        !write(*,*)'Warning Cb (upwind) below zero at node ',i,' with value ',upc(i)
        !		write(100,*)'Warning Cb (upwind) below zero at node ',i,' with value ',upc(i)
		        upc(i) = 0.00001
	        else if(upc(i).gt.0.99999) then
    	        !write(*,*)'Warning Cb (upwind) above one at node ',i,' with value ',upc(i)
    	        !write(100,*)'Warning Cb (upwind) above one at node ',i,' with value ',upc(i)
		        upc(i) = 0.99999
            end if
        
            ! Same for isotopes
	        if(iup(i).lt.0.00001) then
	            !write(*,*)'Warning Ib below zero at node ',i,' with value ',ia(i)
		        !write(100,*)'Warning Ib below zero at node ',i,' with value ',ia(i)
		        iup(i) = 0.00001
            end if	
            
            if(imup(i).lt.0.00001) then
	            !write(*,*)'Warning Imb below zero at node ',i,' with value ',ima(i)
		        !write(100,*)'Warning Imb below zero at node ',i,' with value ',ima(i)
		        imup(i) = 0.00001
            end if
       end do
       
       ccsum=0
       upsum=0
       iupsum=0
       iccsum=0
       uperror=0
       iuperror=0
       
       imcsum=0
       imuperror=0
       imupsum=0
       
       do i = BaseFlagP,TopFlagP
            upsum = upsum + upc(i) 
            ccsum = ccsum + cc(i)
            iupsum = iupsum + iup(i)
            iccsum = iccsum + ic(i)
            imupsum = imupsum + imup(i)
            imcsum = imcsum + imc(i)
       end do
       
       uperror = (upsum+outcomp-(ccsum+outcc))/(ccsum+outcc)
       iuperror = (iupsum+outisot-(iccsum+outci))/(iccsum+outci)
       imuperror= (imupsum+outimt-(imcsum+outcim))/(imcsum+outcim)
       
        do i = BaseFlagP,TopFlagP
           if((pb(i).gt.0))then   
                ca(i) = upc(i)/(1+uperror)
                ia(i) = iup(i)/(1+iuperror)
                ima(i) = imup(i)/(1+iuperror)
           else
               ca(i)=upc(i)
               ia(i)=iup(i)
               ima(i)=imup(i)
           end if
        end do
        cerror = uperror*100
        ierror = iuperror*100
        imerror = imuperror*100
    end if
    
        

    
    
    
    

end subroutine 
!========================================================
!Isotope Solver
!subroutine IsoSolve(z1,z2,dt,dz,pa,pb,ia,ib,ws,wl,is,il,diffdum,maxvel,maxdiff,icc,icerror,z1,z2,advecttype)
!implicit none
!
!integer i,n
!parameter (n=1000000)
!integer z1,z2
!
!real(8) dt,dz,lax
!real(8) pa(-n:n),pb(-n:n),ws(-n:n),wl(-n:n)
!real(8) ia(-n:n),ib(-n:n)
!real(8) il(-n:n),is(-n:n)

!do i = z1+1,z2-2

!	if (i.eq.z2) then
!	ia(i) = ib(i)&
!		&+(dt/dz)*(-il(i)*wl(i)*pb(i) - il(i-1)*ws(i)*(1.0-pb(i)))
!
!	elseif (i.eq.z1) then
!	ia(i) = ib(i)&
!		&-(dt/dz)*(is(i+1)*ws(i+1)*(1.0-pb(i+1)))&
!		&+(dt/dz)*(il(i)*ws(i+1)*(1.0-pb(i+1)))
!	else
! handle all intermediate nodes
! this version uses lax with a very heavy weighting to reduce the diffusion term
! Only update trace element comp if porosity greater than zero 

!if (pb(i).gt.0)then
!	    ia(i) = ((1/(2*dz))*(2*(dz-lax)*ib(i)+lax*ib(i+1)+lax*ib(i-1)))&
!		&-(dt/(2*dz))*(is(i+1)*ws(i+1)*(1.0-pb(i+1)) - is(i-1)*ws(i-1)*(1.0-pb(i-1)))&
!		&+(dt/(2*dz))*(il(i+1)*ws(i+1)*(1.0-pb(i+1)) - il(i-1)*ws(i-1)*(1.0-pb(i-1)))
!	end if
!else
!    ia(i) = ib(i)
!endif
!end do

! this term possibly violates conservation of mass but may be required to stop excess depletion
!	if(ia(i).lt.0.0001) then
!		ia(i) = 0.00001
!	else if(ca(i).gt.1) then
!		ca(i) = 1.0
!	end if
	
!ia(z1) = ib(z1)
!ia(z2) = ib(z2)

!end subroutine 
!========================================================
!TEMPSOLVE
subroutine tempsolve(z1,z2,ha,ta,pa,ca,Hs,Hl,Ts,Tsu, Tl,Lf,Cp)
implicit none

integer i,n,z1,z2
parameter (n=1000000)
real(8) Hs,Hl,Tsu,Tl
real(8) ta(-n:n),ha(-n:n),pa(-n:n),ca(-n:n), ts(-n:n)
real(8) enth_nondim_to_dim,temp_dim_to_nondim,T,H
real(8) Cp,Lf,Solidus,He

do i = z1,z2
	
	H = enth_nondim_to_dim(ha(i),Hs,Hl)
    
    He=Cp*Ts(i)+Lf*ca(i)

	T = (H - Lf*pa(i))/Cp

	if(pa(i).gt.0) then
        if (H.lt.He) then
		    T=Solidus(ca(i),Ts(i))
		end if
	!elseif(pa(i).eq.0) then
	!	T = H/Cp
	end if

	ta(i) = temp_dim_to_nondim(T,Tsu,Tl)

end do
end subroutine
!========================================================
subroutine porossolve(z1,z2,ha,pa,ca,Hs,Hl,Lf,Cp,Tl,Ts,Ae,A1,B1,C1,A2,B2,C2,crysum)

implicit none

 

integer i,k,n,z1,z2

parameter (n=1000000)

real(8) Hs,Hl,Tl,He

real(8) ha(-n:n),pa(-n:n),ca(-n:n), a1(-n:n), b1(-n:n), c1(-n:n), ts(-n:n)

real(8) enth_nondim_to_dim,enth_dim_to_nondim,H, Liquidus,HLiquidus,HSolidus

real(8) Cp,Lf,Ae,A2,B2,C2

real(8) TEST

real(8) sqrt

real(8) fx,fdashx,crysum

 

crysum = 0

 

do i = z1,z2

 

      H = enth_nondim_to_dim(ha(i),Hs,Hl)

! Check if enthalpy is above liquidus: set porosity = 1 if it is

      if(ha(i).gt.enth_dim_to_nondim(HLiquidus(ca(i),Cp,Ae,Lf,A1(i),B1(i),C1(i),A2,B2,C2),Hs,Hl)) then

            TEST = 1

            goto 300

      end if

 

! Check if enthalpy is below solidus: set porosity = 0 if it is

      if(ha(i).lt.enth_dim_to_nondim(HSolidus(ca(i),Ts(i),Cp),Hs,Hl)) then

            TEST = 0

            goto 300

      end if

 

! Check if melting is at solidus (i.e. enthalpy is less than cpT + Lcb).  Calculate porosity as a linear function of enthalpy if it is

! Calculate the enthalpy at the end of eutectic melting for this bulk composition

 

    if (ca(i).lt.Ae) then

	He = Cp*Ts(i) + Lf*ca(i)/Ae

    else

	He = Cp*Ts(i) + Lf*(1-ca(i))/(1-Ae)
    end if

    if(H.le.He)then

        TEST = ca(i)*(H - Hs)/(He - Hs)

!       write(*,*)TEST, ca(i),ha(i),Hl,He,Hs

        goto 300

    endif


! Otherwise solve the cubic function for porosity formed by equating equations (13) and (8)in Solano et al 2014 (express both in terms of T, equate, re-arrange to yield polynomial with 0 on LHS)

 

      if (ca(i).lt.Ae) then

 

            TEST = 1.0

 

! There is an iteration here to update porosity pa based on composition.  Does not use velocity. 

 

600   do k = 1,25

                  fx = (Lf/Cp)*TEST**3 + (C1(i)-H/Cp)*TEST**2 + (B1(i)*ca(i))*TEST + A1(i)*ca(i)*ca(i)

                  fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C1(i)-H/Cp)*TEST + (B1(i)*ca(i))

                  TEST = TEST - fx/fdashx

            end do

 

      else if(ca(i).ge.Ae) then

 

          TEST = 1.0

 

            do k = 1,25

 

            fx = (Lf/Cp)*TEST**3 + (C2-H/Cp + A2 + B2)*TEST**2 + (2*A2*ca(i) - 2*A2 + B2*ca(i) - B2)*TEST + A2*(ca(i)*ca(i)-2*ca(i) + 1)

 

            fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C2-H/Cp + A2 + B2)*TEST + (2*A2*ca(i) - 2*A2 + B2*ca(i) - B2)

 

            !fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C-H/Cp + A + B)*TEST + (B*ca(i)-B)

                                                

                  TEST = TEST - fx/fdashx

 

                  if(TEST.ge.1.0) then

                        TEST = 1.0

                        goto 300

                  end if

 

            end do

 

      end if

 

300   pa(i) = TEST

 

!     if(pa(i).lt.0.00001) then

!           pa(i) = 0

!     end if

crysum = crysum + (1 - pa(i))

 

end do

 

end subroutine

!========================================================
!VELOCITYSOLVE
subroutine velocitysolve(z1,z2,dz,alpha,beta,pa,ws,wl,cla,max_melt_viscosity,min_melt_viscosity,cee,SLT,denscont,drho,aSio2,bSio2,cSio2,dSio2, melt_seg)
implicit none
integer i,n,z1,z2,z0a,z0b
parameter (n=1000000)
real(8) vela,velb,velc,veld,dz,alpha,beta,cee,pad,SLT,velao,velbo,velco
real(8) velbi(-n:n),velci(-n:n),veldi(-n:n),pa(-n:n),ws(-n:n),wl(-n:n),cla(-n:n),denscont(-n:n)
real(8) melt_seg(-n:n)
real(8) max_melt_viscosity, min_melt_viscosity, dimmelvisc, drho,aSio2,bSio2,cSio2,dSio2,ndsio

	do i =  z1,z2
		ws(i) = 0
		wl(i) = 0
	end do


	do i = z1,z2
        if (Melt_seg(i).gt.0) then

        ! Set viscosity function
	        dimmelvisc = 10**(ndsio(cla(i),aSio2,bSio2,cSio2,dSio2)*log10(max_melt_viscosity/min_melt_viscosity))

            if(cee.gt.0)then
        ! Use my equation based on Steve's
        ! Limit maximum porosity used in bulk viscosity equation
                if(pa(i).lt.SLT)then
                    pad=pa(i)
                else
                    pad=SLT
                endif
		        vela = ((pa(i)**alpha)*((1-cee*pad)**beta))/(dz**2) - (cee*beta/(2*(dz**2)))*((pa(i)**alpha)*((1-cee*pad)**(beta-1)))*(pa(i)-pa(i-1))
                velb = -(2/dz**2)*((pa(i)**alpha)*((1-cee*pad)**beta)) - dimmelvisc
		        velc = ((pa(i)**alpha)*((1-cee*pad)**beta))/(dz**2) + (cee*beta/(2*(dz**2)))*((pa(i)**alpha)*((1-cee*pad)**(beta-1)))*(pa(i)-pa(i-1))
            else
        ! Use James' original approach
		        vela = (pa(i)**(alpha+beta))/(dz**2) + (beta/(2*(dz**2)))*(pa(i)**((alpha+beta)-1))*(pa(i)-pa(i-1))
		        velb = -(2/dz**2)*(pa(i)**(alpha+beta)) - dimmelvisc 
		        velc = (pa(i)**(alpha+beta))/(dz**2) - (beta/(2*(dz**2)))*(pa(i)**((alpha+beta)-1))*(pa(i)-pa(i-1))
            endif
        
		    veld = (pa(i)**alpha)*(1-pa(i))*denscont(i)/drho

		    if(i.eq.z2) then
			    velb = velb+velc
		    end if

		    if(i.eq.z1) then
			    velbi(i) = 1.0
			    velci(i) = velc/velb
			    veldi(i) = veld/velb
		    else
			    velbi(i) = 1.0
			    velci(i) = velc/(velb - velci(i-1)*vela)
			    veldi(i) = (veld - veldi(i-1)*vela)/(velb - velci(i-1)*vela)
            end if
            
        end if
	end do

!	ws(z1) = 0
!	wl(z1) = 0

!	ws(z2) = veldi(z2)

	ws(z2+1) = 0
	wl(z2+1) = 0
	
	do i = z2,z1,-1
        if (Melt_seg(i).gt.0) then
		    if (veldi(i).eq.0) then
		        ws(i) = 0 
		    else
		        ws(i) = veldi(i) - velci(i)*ws(i+1)
		    end if
            if (pa(i).gt.0) then
			    wl(i) = ws(i)-ws(i)/pa(i)
            end if
        end if
    end do
	
!	do i = z2,z1,-1
!		if (ws(i).gt.0) then
!			ws(i) = 0
!		endif
!		if (pa(i).gt.0) then
!			wl(i) = ws(i)-ws(i)/pa(i)
!		end if
!	end do

end subroutine
!========================================================
!Melting Rate
subroutine meltrate(z1,z2,dt,dz,pa,pb,ws,melt)
implicit none

integer z1,z2,n,i
parameter (n=1000000)
real(8) dt,dz
real(8) pa(-n:n),pb(-n:n),ws(-n:n)
real(8) melt(-n:n)

do i = z1,z2

		melt(i) = (pa(i) - pb(i))/dt + (pb(i+1)*ws(i+1)-pb(i)*ws(i))/dz

end do

end subroutine

!========================================================
function Solidus(A,Ts)
implicit none
real(8) A,Solidus,Ts

	Solidus = Ts

end function
!========================================================
function Liquidus(Cbulk,Ae,A1,B1,C1,A2,B2,C2)
implicit none
real(8) Ae,Cbulk,Liquidus,A1,B1,C1,A2,B2,C2

if (Cbulk.gt.Ae) then
	Liquidus = A2*(Cbulk**2) + B2*Cbulk + C2
else
	Liquidus = A1*(Cbulk**2) + B1*Cbulk + C1
end if

end function
!========================================================
function SolidComp(T,Cbulk,Ae,A1,B1,C1,A2,B2,C2,Ts)
implicit none
real(8) T,Cbulk,Ae,Solidcomp,Tl,A1,B1,C1,A2,B2,C2,Ts
real(8) Solidus,Liquidus

if (T.gt.Liquidus(Cbulk,Ae,A1,B1,C1,A2,B2,C2)) then
	if (Cbulk.gt.Ae) then
		SolidComp = 1.0
	else
		SolidComp = 0.0	
	end if
else if (T.lt.Solidus(Cbulk,Ts)) then
	SolidComp = Cbulk
else
	if (Cbulk.gt.Ae) then
		SolidComp = 1.0
	else
		SolidComp = 0.0	
	end if
end if

end function
!========================================================
function LiquidComp(T,Cbulk,Ts,Ae,A1,B1,C1,A2,B2,C2)
implicit none
integer i
real(8) T,LC,Cbulk,Ae,Liquidus,Solidus,LiquidComp,Ts,Tl
real(8) A1,B1,C1,A2,B2,C2

	if (T.gt.Liquidus(Cbulk,Ae,A1,B1,C1,A2,B2,C2)) then
		LiquidComp = Cbulk
	else if (T.lt.Solidus(Cbulk,Ts)) then
		LiquidComp = Ae
	else
		if (Cbulk.lt.Ae) then
			LiquidComp = (-B1 - sqrt((B1**2)-4*A1*(C1-T)))/(2*A1)
		else
			LiquidComp = (-B2 + sqrt((B1**2)-4*A2*(C2-T)))/(2*A2)
		end if
	end if

end function
!=========================================================
function HLiquidus(Cbulk,Cp,Ae,Lf,A1,B1,C1,A2,B2,C2)
implicit none
real(8) x,Ae,Cbulk,HLiquidus
real(8) Cp,Lf,A1,B1,C1,A2,B2,C2

!Ae = 0.6417112
!Cp = 1000
!Lf = 550000

if (Cbulk.gt.Ae) then
	HLiquidus = (A2*(Cbulk**2)+B2*Cbulk+C2)*Cp+Lf
else
	HLiquidus = (A1*(Cbulk**2)+B1*Cbulk+C1)*Cp+Lf
end if

end function
!======================================================
function HSolidus(Cbulk,Ts,Cp)
implicit none
real(8) Cbulk,HSolidus,Ts
real(8) Cp,Lf

HSolidus = Ts*Cp

end function

!========================================================
!Functions Converting between Dimensional and Non-Dimensional Temperatures
!========================================================
function temp_dim_to_nondim(T,Ts,Tl)
implicit none
real(8) temp_dim_to_nondim
real(8) T,Ts,Tl

temp_dim_to_nondim = (T - Ts)/(Tl-Ts)

end function
!========================================================
function temp_nondim_to_dim(TPrime,Ts,Tl)
implicit none
real(8) temp_nondim_to_dim
real(8) TPrime,Ts,Tl

temp_nondim_to_dim = TPrime*(Tl-Ts) + Ts

end function
!========================================================
function enth_dim_to_nondim(H,Hs,Hl)
implicit none
real(8) enth_dim_to_nondim
real(8) H,Hs,Hl

enth_dim_to_nondim = (H - Hs)/(Hl-Hs)

end function
!========================================================
function enth_nondim_to_dim(HPrime,Hs,Hl)
implicit none
real(8) enth_nondim_to_dim
real(8) HPrime,Hs,Hl

enth_nondim_to_dim = HPrime*(Hl-Hs) + Hs

end function
!========================================================
!Functions Describing the equilibrium isotope concentration
!========================================================
function IsoSolid(Phi,Conc,BulkComp,LiqComp,SolComp,PartitionA,PartitionB)
implicit none
real(8) IsoSolid
real(8) conc,PartitionA,PartitionB,Phi
real(8) BulkPartition
real(8) BulkComp,LiqComp,SolComp
real(8) PPP,DDD, Ps, Pl

!PPP = PartitionA*(1-LiqComp)+ PartitionB*LiqComp
!DDD = (1-SolComp)*PartitionA + SolComp*PartitionB

Ps = (1-phi)*((1-SolComp)*PartitionA + SolComp*PartitionB)
Pl = phi*(PartitionA*(1-LiqComp)+ PartitionB*LiqComp)

if(phi.eq.1.0) then
	IsoSolid = 0.0
else
!	IsoSolid = Conc*(DDD-PPP*Phi)/((1-Phi)*(DDD+Phi*(1-PPP)))
    IsoSolid = Conc*(Ps-phi*Pl)/((1-phi)*(Ps+phi*(1-Pl)))
end if

if(IsoSolid.lt.0.0) then
	Isosolid = 0.0
end if

end function
!========================================================
function IsoLiquid(Phi,Conc,BulkComp,LiqComp,SolComp,PartitionA,PartitionB)
implicit none
real(8) IsoLiquid
real(8) conc,PartitionA,PartitionB,Phi
real(8) BulkPartition
real(8) BulkComp,LiqComp,SolComp
real(8) PPP,DDD

PPP = PartitionA*(1-LiqComp)+ PartitionB*LiqComp
DDD = (1-SolComp)*PartitionA + SolComp*PartitionB

if(Phi.eq.0.0) then
	IsoLiquid = 0.0
else
	IsoLiquid = Conc/(DDD + Phi*(1-PPP))
end if

if(IsoLiquid.lt.0.0) then
	IsoLiquid = 0.0
end if

end function
!========================================================
subroutine RandomNumber(Random,depthnodes,Seed,empexp)
implicit none

integer random
real i,j                 ! Counts random numbers
real(8) ran
integer depthnodes,empexp
integer seed

j = ran(seed)
i = 2*j-1

random = int(-depthnodes*(i**empexp)/2)
!write(*,*)j,i, random, depthnodes
!write(100,*) 'RANDOM',j,i, random, depthnodes, empexp

end subroutine

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine lowershiftmass(j1,j2,z1,z2,oldj3,oldj4,ovtn,TransNodes,radg,DMFn,ha,ta,&
                    &pa,ca,cd,csa,cla,ica,icc,icsa,icla,ima,imc,imsa,imla,hb,tb,&
                    &pb,cb,csb,clb,icb,icsb,iclb,imb,imsb,imlb,clostim,density,&
                    & PartitionA,PartitionB,Avca,Avcd,AvT,Avh,Avcl,Avcs,Avpa,Avdensill,&
                    &runtime,year,tau,baseflagP,topflagP,dz,delta,asio2,bsio2,csio2,dsio2,&
                    &AvA1, AvB1, AvC1, AvTs,A1,B1,C1,Ts,AvIa,AvIma,AvIcc,AvImc,Tl,Tsu,Closetemp,midD,OGMIDLIM,OGUPLIM,d21,melt_seg,pmag)
implicit none
integer i,n,j1,j2,z1,z2,DMFn,baseflagP,topflagP,TransNodes,STransNodes, ovtn
integer d21
real(8) radg, dz, delta,asio2,bsio2,csio2,dsio2
parameter (n=1000000)


real(8) ha(-n:n)	    !Enthalpy
real(8) ta(-n:n)	    !Temperature
real(8) pa(-n:n)	    !Porosity
real(8) ca(-n:n)	    !Bulk Composition
real(8) cd(-n:n)        !Control composition
real(8) csa(-n:n)	    !Solid Composition
real(8) cla(-n:n)	    !Liquid Composition
real(8) ica(-n:n)	    !Isotope Bulk Composition
real(8) icc(-n:n)
real(8) icsa(-n:n)	    !Solid Composition
real(8) icla(-n:n)	    !Liquid Composition
real(8) ima(-n:n)	    !Isotope Bulk Composition
real(8) imsa(-n:n)	    !Solid Composition
real(8) imla(-n:n)	    !Liquid Composition
real(8) clostim(-n:n)	!Closure time
real(8) density(-n:n)   !Density
real(8) hb(-n:n)	    !Enthalpy
real(8) tb(-n:n)	    !Temperature
real(8) pb(-n:n)	    !Porosity
real(8) cb(-n:n)	    !Bulk Composition
real(8) csb(-n:n)	    !Solid Composition
real(8) clb(-n:n)	    !Liquid Composition
real(8) icb(-n:n)	    !Isotope Bulk Composition
real(8) icsb(-n:n)	    !Solid Composition
real(8) iclb(-n:n)	    !Liquid Composition
real(8) imb(-n:n)	    !Isotope Bulk Composition
real(8) imsb(-n:n)	    !Solid Composition
real(8) imlb(-n:n)	    !Liquid Composition
real(8) imc(-n:n)
real(8) A1(-n:n)
real(8) B1(-n:n)
real(8) C1(-n:n)
real(8) Ts(-n:n)
real(8) pmag(-n:n)

real(8) melt_seg(-n:n)
real(8) melt_seg1(-n:n)



real(8) ha1(-n:n)	    !Enthalpy
real(8) ta1(-n:n)	    !Temperature
real(8) pa1(-n:n)	    !Porosity
real(8) ca1(-n:n)	    !Bulk Composition
real(8) cd1(-n:n)       !Control composition
real(8) csa1(-n:n)	    !Solid Composition
real(8) cla1(-n:n)	    !Liquid Composition
real(8) ica1(-n:n)	    !Isotope Bulk Composition
real(8) icc1(-n:n)
real(8) icsa1(-n:n)	    !Solid Composition
real(8) icla1(-n:n)	    !Liquid Composition
real(8) ima1(-n:n)	    !Isotope Bulk Composition
real(8) imc1(-n:n)     
real(8) imsa1(-n:n)	    !Solid Composition
real(8) imla1(-n:n)	    !Liquid Composition
real(8) clostim1(-n:n)	!Closure time
real(8) density1(-n:n)

real(8) hb1(-n:n)	    !Enthalpy
real(8) tb1(-n:n)	    !Temperature
real(8) pb1(-n:n)	    !Porosity
real(8) cb1(-n:n)	    !Bulk Composition
real(8) csb1(-n:n)	    !Solid Composition
real(8) clb1(-n:n)	    !Liquid Composition
real(8) icb1(-n:n)	    !Isotope Bulk Composition
real(8) icsb1(-n:n)	    !Solid Composition
real(8) iclb1(-n:n)	    !Liquid Composition
real(8) imb1(-n:n)	    !Isotope Bulk Composition
real(8) imsb1(-n:n)	    !Solid Composition
real(8) imlb1(-n:n)	    !Liquid Composition
real(8) A11(-n:n)
real(8) B11(-n:n)
real(8) C11(-n:n)
real(8) Ts1(-n:n)
real(8) pmag1(-n:n)

real(8) PartitionA,PartitionB

real(8) SolidComp,LiquidComp
real(8) IsoSolid,IsoLiquid
real(8) temp_nondim_to_dim

real(8)Avcd,Avca,AvT,Avh,Avcl,Avcs,Avpa,Avdensill,Avic,icasum,imasum,AvIm
real(8)runtime,year,tau

real(8) AvA1, AvB1, AvC1, AvTs, CloseTemp,Tl
real(8) AvIcc, avia, tsu, AvImc, AvIma

integer midD, oldj3, oldj4
integer ogmidlim, oguplim

STransNodes=int((TransNodes)*(radg))

!! Saving values

        do i=z1,z2
            ha1(i) = ha(i)
            ta1(i) = ta(i)
            pa1(i) = pa(i)
		
            ca1(i) = ca(i)
            cd1(i) = cd(i)
            csa1(i) = csa(i)
            cla1(i) = cla(i)
		
            ica1(i) = ica(i)
            icc1(i) = icc(i)
            icsa1(i) = icsa(i)
            icla1(i) = icla(i)

            ima1(i) = ima(i)
            imc1(i) = imc(i)
            imsa1(i) = imsa(i)
            imla1(i) = imla(i)
		
            clostim1(i) = clostim(i)
            density1(i) = density(i)
            
            hb1(i) = hb(i)
            tb1(i) = tb(i)
            pb1(i) = pb(i)
		
            cb1(i) = cb(i)
            csb1(i) = csb(i)
            clb1(i) = clb(i)
		
            icb1(i) = icb(i)
            icsb1(i) = icsb(i)
            iclb1(i) = iclb(i)

            imb1(i) = imb(i)
            imsb1(i) = imsb(i)
            imlb1(i) = imlb(i)
            
            A11(i)=A1(i)
            B11(i)=B1(i)
            C11(i)=C1(i)
            Ts1(i)=Ts(i)
            pmag1(i) = pmag(i)
            
            Melt_seg1(i)=Melt_seg(i)
           
        end do

  
!! Shifting - below LCMR
        do i=z1,j1-Transnodes
 		
		    ha(i-(STransNodes-Transnodes)) = ha1(i)
            
		    ta(i-(STransNodes -Transnodes)) = ta1(i)
            
		    pa(i-(STransNodes -Transnodes)) = pa1(i)
		
		    ca(i-(STransNodes -Transnodes)) = ca1(i)
            cd(i-(STransNodes -Transnodes))= cd1(i)
		    csa(i-(STransNodes -Transnodes)) = csa1(i)
		    cla(i-(STransNodes -Transnodes)) = cla1(i)
		
		    ica(i-(STransNodes -Transnodes))= ica1(i)
            icc(i-(STransNodes -Transnodes)) = icc(i)
		    icsa(i-(STransNodes -Transnodes))= icsa1(i)
		    icla(i-(STransNodes -Transnodes))= icla1(i)

		    ima(i-(STransNodes -Transnodes)) = ima1(i)
            imc(i-(STransNodes -Transnodes)) = imc1(i)
		    imsa(i-(STransNodes -Transnodes)) = imsa1(i)
		    imla(i-(STransNodes -Transnodes)) = imla1(i)
            
            hb(i-(STransNodes -Transnodes))= hb1(i)
		    tb(i-(STransNodes -Transnodes)) = tb1(i)
		    pb(i-(STransNodes -Transnodes)) = pb1(i)
		
		    cb(i-(STransNodes -Transnodes)) = cb1(i)
		    csb(i-(STransNodes -Transnodes)) = csb1(i)
		    clb(i-(STransNodes -Transnodes)) = clb1(i)
		
		    icb(i-(STransNodes -Transnodes)) = icb1(i)
		    icsb(i-(STransNodes -Transnodes))= icsb1(i)
		    iclb(i-(STransNodes -Transnodes)) = iclb1(i)

		    imb(i-(STransNodes -Transnodes)) = imb1(i)
		    imsb(i-(STransNodes -Transnodes)) = imsb1(i)
		    imlb(i-(STransNodes -Transnodes)) = imlb1(i)
		    
            density(i-(STransNodes -Transnodes)) = density1(i)
            pmag(i-(STransNodes -Transnodes)) = pmag1(i)
            
            clostim(i-(STransNodes -Transnodes))= clostim1(i)
            
            A1(i-(STransNodes -Transnodes))= A11(i)
            B1(i-(STransNodes -Transnodes))= B11(i)
            C1(i-(STransNodes -Transnodes))= C11(i)
            Ts(i-(STransNodes -Transnodes))= Ts1(i)
           
            Melt_seg(i-(STransNodes -Transnodes))= Melt_seg1(i)
        end do
        
        z1=z1-(STransNodes -Transnodes)
        
        if (d21.le.(j1-Transnodes)) then 
            d21=d21-(STransNodes -Transnodes)
        else if ((d21.gt.(j1-Transnodes)).AND.((d21.le.(DMFn)))) then 
            d21=d21-STransNodes
        end if
       
        
        if (OGmidLIM.le.(j1-Transnodes)) then 
            OGmidLIM=OGmidLIM-(STransNodes -Transnodes)
        else if ((OGmidLIM.gt.(j1-Transnodes)).AND.((OGmidLIM.le.(DMFn)))) then 
            OGmidLIM=OGmidLIM-STransNodes
        end if
        
        if (OGupLIM.le.(j1-Transnodes)) then 
            OGupLIM=OGupLIM-(STransNodes -Transnodes)
        else if ((OGupLIM.gt.(j1-Transnodes)).AND.((OGupLIM.le.(DMFn)))) then 
            OGupLIM=OGupLIM-STransNodes
        end if
        
        if (((DMFn-STransnodes).le.midD).AND.((DMFn-Stransnodes).gt.OGmidLim)) then
            midD=DMfn-STransNodes
        end if
        
       ! baseflagP=baseflagP-(STransNodes -Transnodes)
    ! write(*,*) midD
        !topflagP=topflagP-(STransNodes -Transnodes)
!! Shifitng - from intrusion site down to LCMR
        
        do i=j1+1, DMFn
            
            ha(i-STransNodes ) = ha1(i)
		    ta(i-STransNodes ) = ta1(i)
		    pa(i-STransNodes ) = pa1(i)
		
		    ca(i-STransNodes ) = ca1(i)
            cd(i-STransNodes ) = cd1(i)
		    csa(i-STransNodes ) = csa1(i)
		    cla(i-STransNodes ) = cla1(i)
		
		    ica(i-STransNodes ) = ica1(i)
            icc(i-STransNodes ) = icc1(i)
		    icsa(i-STransNodes ) = icsa1(i)
		    icla(i-STransNodes ) = icla1(i)

		    ima(i-STransNodes ) = ima1(i)
            imc(i-STransNodes) = imc1(i)
		    imsa(i-STransNodes ) = imsa1(i)
		    imla(i-STransNodes ) = imla1(i)
            
            hb(i-STransNodes ) = hb1(i)
		    tb(i-STransNodes ) = tb1(i)
		    pb(i-STransNodes ) = pb1(i)
		
		    cb(i-STransNodes ) = cb1(i)
		    csb(i-STransNodes ) = csb1(i)
		    clb(i-STransNodes ) = clb1(i)
		
		    icb(i-STransNodes ) = icb1(i)
		    icsb(i-STransNodes ) = icsb1(i)
		    iclb(i-STransNodes ) = iclb1(i)

		    imb(i-STransNodes ) = imb1(i)
		    imsb(i-STransNodes ) = imsb1(i)
		    imlb(i-STransNodes ) = imlb1(i)
            pmag(i-STransNodes ) = pmag1(i)
		
            clostim(i-STransNodes ) = clostim1(i)
        
            density(i-STransNodes )= density1(i)
            
            A1(i-STransNodes )= A11(i)
            B1(i-STransNodes )= B11(i)
            C1(i-STransNodes )= C11(i)
            TS(i-STransNodes )= TS1(i)
            
            Melt_seg(i-STransNodes )= Melt_seg1(i)
            
        end do

        if ((oldj3.ge.(j1+1)).AND.((oldj3.le.(DMFn)))) then
            oldj3=oldj3-STransNodes
        end if 
        
        if ((oldj4.ge.(j1+1)).AND.((oldj4.le.(DMFn)))) then
            oldj4=oldj4-STransNodes
        end if 
!!! INTUDE THE MATERIAL
     
        do i=DMFn-(STransNodes-1), DMFn
            ha(i) = Avh
            ta(i) = AvT
            pa(i) = AvPa
            ca(i) = AvCa
            !cd(i) = Avcd
            
            
            csa(i) = Avcs
            cla(i) = Avcl
            
            ica(i) = AvIa
            icsa(i) = IsoSolid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            icla(i) = IsoLiquid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            
            ima(i) = AvIma
            imsa(i) = IsoSolid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            imla(i) = IsoLiquid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            
            A1(i)=AvA1
            B1(i)=AvB1
            C1(i)=AvC1
            Ts(i)=AvTs
            
            Melt_seg(i)=1
            pmag(i) = 1

            if(((AvT*(Tl-Tsu)+Tsu).ge.CloseTemp))then
                clostim(i) = runtime*tau/year/1E6
            else
                clostim(i)=-999
            end if     
        end do
        
        do i= DMFn-(Transnodes-1), DMFn
            cd(i)=Avcd
            icc(i)=Avicc
            imc(i)=Avimc
        end do
        do i=DMFn-(STransNodes -1), DMFn-(Transnodes-1)-1
            cd(i)=Avca
            icc(i)=Avia
            imc(i)=Avima
        end do
        



!j1=j1-Transnodes-(STransNodes-Transnodes)

!j2=j2-Transnodes-(STransNodes-Transnodes)










end subroutine    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine eruptandshiftmass(j1b,j2b,z1,z2,oldL3,oldL4,DMFn,AmountI,AILower,AELower,Avh&
                    &, AvT, AvPa, AvCa, Avcs, Avcl,ha,ta,pa,ca,cd,csa,&
                    &cla,ica,icc,icsa,icla,ima,imc,imsa,imla,hb,tb,pb,cb,csb,&
                    &clb,icb,icsb,iclb,imb,imsb,imlb,clostim,density,&
                    & PartitionA,PartitionB,&
                    &runtime,year,tau,outcc,outcomp,baseflagP,topflagP,AvA1,&
                    &AvB1, AvC1, AvTs,A1,B1,C1,Ts,&
                    &AvIa, AvIma, outci,outcim, outisot,outimt,Tl,tsu,&
                    &CloseTemp,upperD,midD,OGMIDLIM, OGUPLIM,d21,dmfnL,melt_seg_flag,&
                    &melt_seg, pmag)
implicit none
integer i,n,j1b,j2b,z1,z2,baseflagP,topflagP,AmountI, AILower, AELower,DMFn
integer oldL3, oldL4
real(8) outcc, outcomp, outcc1,outcomp1
parameter (n=1000000)


real(8) ha(-n:n)	    !Enthalpy
real(8) ta(-n:n)	    !Temperature
real(8) pa(-n:n)	    !Porosity
real(8) ca(-n:n)	    !Bulk Composition
real(8) cd(-n:n)        !Control composition
real(8) csa(-n:n)	    !Solid Composition
real(8) cla(-n:n)	    !Liquid Composition
real(8) ica(-n:n)	    !Isotope Bulk Composition
real(8) icc(-n:n)
real(8) icsa(-n:n)	    !Solid Composition
real(8) icla(-n:n)	    !Liquid Composition
real(8) ima(-n:n)	    !Isotope Bulk Composition
real(8) imc(-n:n)
real(8) imsa(-n:n)	    !Solid Composition
real(8) imla(-n:n)	    !Liquid Composition
real(8) clostim(-n:n)	!Closure time
real(8) density(-n:n)   !Density
real(8) hb(-n:n)	    !Enthalpy
real(8) tb(-n:n)	    !Temperature
real(8) pb(-n:n)	    !Porosity
real(8) cb(-n:n)	    !Bulk Composition
real(8) csb(-n:n)	    !Solid Composition
real(8) clb(-n:n)	    !Liquid Composition
real(8) icb(-n:n)	    !Isotope Bulk Composition
real(8) icsb(-n:n)	    !Solid Composition
real(8) iclb(-n:n)	    !Liquid Composition
real(8) imb(-n:n)	    !Isotope Bulk Composition
real(8) imsb(-n:n)	    !Solid Composition
real(8) imlb(-n:n)	    !Liquid Composition
real(8) A1(-n:n)
real(8) B1(-n:n)
real(8) C1(-n:n)
real(8) Ts(-n:n)
real(8) pmag(-n:n)



real(8) ha1(-n:n)	    !Enthalpy
real(8) ta1(-n:n)	    !Temperature
real(8) pa1(-n:n)	    !Porosity
real(8) ca1(-n:n)	    !Bulk Composition
real(8) cd1(-n:n)       !Control composition
real(8) csa1(-n:n)	    !Solid Composition
real(8) cla1(-n:n)	    !Liquid Composition
real(8) ica1(-n:n)	    !Isotope Bulk Composition
real(8) icc1(-n:n)
real(8) icsa1(-n:n)	    !Solid Composition
real(8) icla1(-n:n)	    !Liquid Composition
real(8) ima1(-n:n)	    !Isotope Bulk Composition
real(8) imsa1(-n:n)	    !Solid Composition
real(8) imla1(-n:n)	    !Liquid Composition
real(8) imc1(-n:n)
real(8) clostim1(-n:n)	!Closure time
real(8) density1(-n:n)

real(8) hb1(-n:n)	    !Enthalpy
real(8) tb1(-n:n)	    !Temperature
real(8) pb1(-n:n)	    !Porosity
real(8) cb1(-n:n)	    !Bulk Composition
real(8) csb1(-n:n)	    !Solid Composition
real(8) clb1(-n:n)	    !Liquid Composition
real(8) icb1(-n:n)	    !Isotope Bulk Composition
real(8) icsb1(-n:n)	    !Solid Composition
real(8) iclb1(-n:n)	    !Liquid Composition
real(8) imb1(-n:n)	    !Isotope Bulk Composition
real(8) imsb1(-n:n)	    !Solid Composition
real(8) imlb1(-n:n)	    !Liquid Composition
real(8) A11(-n:n)
real(8) B11(-n:n)
real(8) C11(-n:n)
real(8) Ts1(-n:n)
real(8) pmag1(-n:n)

real(8) melt_seg(-n:n)
real(8) melt_seg1(-n:n)

real(8) PartitionA,PartitionB

real(8) SolidComp,LiquidComp
real(8) IsoSolid,IsoLiquid
real(8) temp_nondim_to_dim

real(8)Avcd,Avca,AvT,Avh,Avcl,Avcs,Avpa,icasum,imasum,AvIma
real(8)runtime,year,tau
real(8) AvA1, AvB1, AvC1, AvTs

real(8) oi, avia, outci1, avicc, outci, outisot
real(8) outcim1, avimc, outcim, outimt
real(8) tl, tsu, CloseTemp
integer upperD, midD
integer ogmidlim, oguplim, d21
integer DMFNL

integer melt_seg_flag

!! Saving values

        do i=z1,z2
            ha1(i) = ha(i)
            ta1(i) = ta(i)
            pa1(i) = pa(i)
		
            ca1(i) = ca(i)
            cd1(i) = cd(i)
            csa1(i) = csa(i)
            cla1(i) = cla(i)
		
            ica1(i) = ica(i)
            icc1(i) = icc(i)
            icsa1(i) = icsa(i)
            icla1(i) = icla(i)

            ima1(i) = ima(i)
            imc1(i) = imc(i)
            imsa1(i) = imsa(i)
            imla1(i) = imla(i)
		
            clostim1(i) = clostim(i)
            density1(i) = density(i)
            
            hb1(i) = hb(i)
            tb1(i) = tb(i)
            pb1(i) = pb(i)
		
            cb1(i) = cb(i)
            csb1(i) = csb(i)
            clb1(i) = clb(i)
		
            icb1(i) = icb(i)
            icsb1(i) = icsb(i)
            iclb1(i) = iclb(i)

            imb1(i) = imb(i)
            imsb1(i) = imsb(i)
            imlb1(i) = imlb(i)
            
            A11(i)=A1(i)
            B11(i)=B1(i)
            C11(i)=C1(i)
            Ts1(i)=Ts(i)
            
            
            Melt_seg1(i)=Melt_seg(i)
            pmag1(i) = pmag(i)
            
            

        end do
        
        outcc1=0.0
        outcomp1=0.0
        outci1=0.0
        outcim1=0.0
        
        do i=j2b,j1b
            outcc1=outcc1+cd(i)
            outci1=outci1+icc(i)
            outcim1=outcim1+imc(i)
            !outcomp1=outcomp1+ca(i)
        end do
        
        Avcd=outcc1/(abs(j1b-j2b)+1)
        Avicc=outci1/(abs(j1b-j2b)+1)
        Avimc=outcim1/(abs(j1b-j2b)+1)
        
        outcc=outcc+Avcd*AELower
        outci=outci+Avicc*AELower
        outcim=outcim+Avimc*AELower
        
        outcomp=outcomp+Avca*AELower
        outisot=outisot+AvIa*AELower
        outimt=outimt+AvIma*AELower
        
!! Shifting - from below the high melt fraction to the bottom of the crust
        
        do i=z1,j2b-1
 		
		    ha(i+(abs(j1b-j2b)+1)-AmountI) = ha1(i)
		    ta(i+(abs(j1b-j2b)+1)-AmountI) = ta1(i)
		    pa(i+(abs(j1b-j2b)+1)-AmountI) = pa1(i)
		
		    ca(i+(abs(j1b-j2b)+1)-AmountI) = ca1(i)
            cd(i+(abs(j1b-j2b)+1)-AmountI)= cd1(i)
		    csa(i+(abs(j1b-j2b)+1)-AmountI) = csa1(i)
		    cla(i+(abs(j1b-j2b)+1)-AmountI) = cla1(i)
		
		    ica(i+(abs(j1b-j2b)+1)-AmountI)= ica1(i)
            icc(i+(abs(j1b-j2b)+1)-AmountI) = icc1(i)
		    icsa(i+(abs(j1b-j2b)+1)-AmountI)= icsa1(i)
		    icla(i+(abs(j1b-j2b)+1)-AmountI)= icla1(i)

		    ima(i+(abs(j1b-j2b)+1)-AmountI) = ima1(i)
            imc(i+(abs(j1b-j2b)+1)-AmountI) = imc1(i)
		    imsa(i+(abs(j1b-j2b)+1)-AmountI) = imsa1(i)
		    imla(i+(abs(j1b-j2b)+1)-AmountI) = imla1(i)
            
            hb(i+(abs(j1b-j2b)+1)-AmountI)= hb1(i)
		    tb(i+(abs(j1b-j2b)+1)-AmountI) = tb1(i)
		    pb(i+(abs(j1b-j2b)+1)-AmountI) = pb1(i)
		
		    cb(i+(abs(j1b-j2b)+1)-AmountI) = cb1(i)
		    csb(i+(abs(j1b-j2b)+1)-AmountI)= csb1(i)
		    clb(i+(abs(j1b-j2b)+1)-AmountI) = clb1(i)
		
		    icb(i+(abs(j1b-j2b)+1)-AmountI) = icb1(i)
		    icsb(i+(abs(j1b-j2b)+1)-AmountI)= icsb1(i)
		    iclb(i+(abs(j1b-j2b)+1)-AmountI) = iclb1(i)

		    imb(i+(abs(j1b-j2b)+1)-AmountI) = imb1(i)
		    imsb(i+(abs(j1b-j2b)+1)-AmountI) = imsb1(i)
		    imlb(i+(abs(j1b-j2b)+1)-AmountI) = imlb1(i)
		    
            density(i+(abs(j1b-j2b)+1)-AmountI) = density1(i)
            
            clostim(i+(abs(j1b-j2b)+1)-AmountI)= clostim1(i)
            
            
            A1(i+(abs(j1b-j2b)+1)-AmountI) = A11(i)
            B1(i+(abs(j1b-j2b)+1)-AmountI) = B11(i)
            C1(i+(abs(j1b-j2b)+1)-AmountI) = C11(i)
            TS(i+(abs(j1b-j2b)+1)-AmountI) = TS1(i)
            
           
           Melt_seg(i+(abs(j1b-j2b)+1)-AmountI) = Melt_seg1(i)
           pmag(i+(abs(j1b-j2b)+1)-AmountI) = pmag1(i)
           
        end do
        
        z1=z1+(abs(j1b-j2b)+1)-AmountI
        
        
        
        midD=midD+(abs(j1b-j2b)+1)-AmountI
        
        
        
        if (d21.le.(j2b-1)) then 
            d21=d21+(abs(j1b-j2b)+1)-AmountI
        else if ((d21.gt.(j2b-1)).AND.((d21.le.(DMFn)))) then 
            d21=d21-AmountI
        end if
        
        if (dMFNL.le.(j2b-1)) then 
            dMFNL=dMFNL+(abs(j1b-j2b)+1)-AmountI
        else if ((dMFNL.gt.(j2b-1)).AND.((dMFNL.le.(DMFn)))) then 
            dMFNL=dMFNL-AmountI
        end if
        
        
        if (OGmidLIM.le.(j2b-1)) then 
            OGmidLIM=OGmidLIM+(abs(j1b-j2b)+1)-AmountI
        else if ((OGmidLIM.gt.(j2b-1)).AND.((OGmidLIM.le.(DMFn)))) then 
            OGmidLIM=OGmidLIM-AmountI
        end if
        
        if (OGupLIM.le.(j2b-1)) then 
            OGupLIM=OGupLIM+(abs(j1b-j2b)+1)-AmountI
        else if ((OGupLIM.gt.(j2b-1)).AND.((OGupLIM.le.(DMFn)))) then 
            OGupLIM=OGupLIM-AmountI
        end if
        
        if (oldL3.ne.5000000) then
            oldL3=oldL3+(abs(j1b-j2b)+1)-AmountI
        end if
        
        if (oldL4.ne.5000000) then
            oldL4=oldL4+(abs(j1b-j2b)+1)-AmountI
        end if
       ! baseflagP=baseflagP+(abs(j1b-j2b)+1)-AmountI
        if (((DMFN-AmountI).le.upperD).AND.((DMFN-AmountI).gt.OGupLIM)) then
            upperD=DMfn-AmountI
        end if
        
        !! Shfiting from intrusion site down to the where the high melt fraction layer once was

        do i=j1b+1, DMFn-1
            
            ha(i-AmountI) = ha1(i)
		    ta(i-AmountI) = ta1(i)
		    pa(i-AmountI) = pa1(i)
		
		    ca(i-AmountI) = ca1(i)
            cd(i-AmountI) = cd1(i)
		    csa(i-AmountI) = csa1(i)
		    cla(i-AmountI) = cla1(i)
		
		    ica(i-AmountI) = ica1(i)
            icc(i-AmountI) = icc1(i)
		    icsa(i-AmountI) = icsa1(i)
		    icla(i-AmountI) = icla1(i)

		    ima(i-AmountI) = ima1(i)
            imc(i-AmountI) = imc1(i)
		    imsa(i-AmountI) = imsa1(i)
		    imla(i-AmountI) = imla1(i)
            
            hb(i-AmountI) = hb1(i)
		    tb(i-AmountI) = tb1(i)
		    pb(i-AmountI) = pb1(i)
		
		    cb(i-AmountI) = cb1(i)
		    csb(i-AmountI) = csb1(i)
		    clb(i-AmountI) = clb1(i)
		
		    icb(i-AmountI) = icb1(i)
		    icsb(i-AmountI) = icsb1(i)
		    iclb(i-AmountI) = iclb1(i)

		    imb(i-AmountI) = imb1(i)
		    imsb(i-AmountI) = imsb1(i)
		    imlb(i-AmountI) = imlb1(i)
            
            A1(i-AmountI) = A11(i)
            B1(i-AmountI) = B11(i)
            C1(i-AmountI) = C11(i)
            Ts(i-AmountI) = Ts1(i)
            
            clostim(i-AmountI) = clostim1(i)
        
            density(i-AmountI)= density1(i)
            
            
            Melt_seg(i-AmountI)= Melt_seg1(i)
            pmag(i-AmountI)=pmag1(i)
            
            
        end do

!!! INTUDE THE MATERIAL
 if (AmountI.gt.0) then       
        do i=DMFn-AmountI, DMFn
            ha(i) = Avh
            ta(i) = AvT
            pa(i) = AvPa
            ca(i) = AvCa
            !cd(i) = Avcd
            
            
            csa(i) = Avcs
            cla(i) = Avcl
		
            ica(i) = AvIa
            icsa(i) = IsoSolid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            icla(i) = IsoLiquid(pa(i),ica(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            !icc(i) = AvIC

            ima(i) = Avima
            imsa(i) = IsoSolid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
            imla(i) = IsoLiquid(pa(i),ima(i),ca(i),cla(i),csa(i),PartitionA,PartitionB)
        
            A1(i) = AvA1
            B1(i) = AvB1
            C1(i) = AvC1
            Ts(i) = AvTs
            pmag(i) = 1
        
		    
            if(((AvT*(Tl-Tsu)+Tsu).ge.CloseTemp))then
                clostim(i) = runtime*tau/year/1E6
            else
                clostim(i)=-999
            end if
        
            if (Melt_seg_flag.gt.0) then
		        Melt_seg(i)=0
            else
                Melt_seg(i)=1
            end if
                        
               
        end do
        
        do i= DMFn-AmountI, DMFn-AmountI+AILower
            cd(i)=Avcd
            icc(i)=Avicc
            imc(i)=Avimc
        end do
        do i=DMFn-AmountI+AILower+1, DMFn
            cd(i)=Avca
            icc(i)=AvIa
            imc(i)=AvIma
        end do
        
end if 







end subroutine
    
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

subroutine porosshift(ha,pa,ca,Hs,Hl,Lf,Cp,Tl,Ts,Ae,A1,B1,C1,A2,B2,C2)

implicit none

 


integer k

real(8) Hs,Hl,Ts,Tl,He

real(8) ha,pa,ca

real(8) enth_nondim_to_dim,enth_dim_to_nondim,H, Liquidus,HLiquidus,HSolidus

real(8) Cp,Lf,Ae,A1,B1,C1,A2,B2,C2

real(8) TEST

real(8) sqrt

real(8) fx,fdashx

 



 


 

      H = enth_nondim_to_dim(ha,Hs,Hl)

! Check if enthalpy is above liquidus: set porosity = 1 if it is

      if(ha.gt.enth_dim_to_nondim(HLiquidus(ca,Cp,Ae,Lf,A1,B1,C1,A2,B2,C2),Hs,Hl)) then

            TEST = 1

            goto 300

      end if

 

! Check if enthalpy is below solidus: set porosity = 0 if it is

      if(ha.lt.enth_dim_to_nondim(HSolidus(ca,Ts,Cp),Hs,Hl)) then

            TEST = 0

            goto 300

      end if

 

! Check if melting is at solidus (i.e. enthalpy is less than cpT + Lcb).  Calculate porosity as a linear function of enthalpy if it is

! Calculate the enthalpy at the end of eutectic melting for this bulk composition

 

    if (ca.lt.Ae) then

	He = Cp*Ts + Lf*ca/Ae

    else

	He = Cp*Ts + Lf*(1-ca)/(1-Ae)
    end if

    if(H.le.He)then

        TEST = ca*(H - Hs)/(He - Hs)

!       write(*,*)TEST, ca(i),ha(i),Hl,He,Hs

        goto 300

    endif


! Otherwise solve the cubic function for porosity formed by equating equations (13) and (8)in Solano et al 2014 (express both in terms of T, equate, re-arrange to yield polynomial with 0 on LHS)

 

      if (ca.lt.Ae) then

 

            TEST = 1.0

 

! There is an iteration here to update porosity pa based on composition.  Does not use velocity. 

 

600   do k = 1,25

                  fx = (Lf/Cp)*TEST**3 + (C1-H/Cp)*TEST**2 + (B1*ca)*TEST + A1*ca*ca

                  fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C1-H/Cp)*TEST + (B1*ca)

                  TEST = TEST - fx/fdashx

            end do

 

      else if(ca.ge.Ae) then

 

          TEST = 1.0

 

            do k = 1,25

 

            fx = (Lf/Cp)*TEST**3 + (C2-H/Cp + A2 + B2)*TEST**2 + (2*A2*ca - 2*A2 + B2*ca - B2)*TEST + A2*(ca*ca-2*ca + 1)

 

            fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C2-H/Cp + A2 + B2)*TEST + (2*A2*ca - 2*A2 + B2*ca - B2)

 

            !fdashx = 3*(Lf/Cp)*TEST**2 + 2*(C-H/Cp + A + B)*TEST + (B*ca(i)-B)

                                                

                  TEST = TEST - fx/fdashx

 

                  if(TEST.ge.1.0) then

                        TEST = 1.0

                        goto 300

                  end if

 

            end do

 

      end if

 

300   pa = TEST

 

!     if(pa(i).lt.0.00001) then

!           pa(i) = 0+1

!     end if



 



 

    end subroutine
    

    
    

function ndsio(ca,aSio2,bSio2,cSio2,dSio2)



implicit none


real(8) ca,ndsio
real(8) aSio2,bSio2,cSio2,dSio2
real(8) MaxSio2, MinSio2,Sio2

Sio2=aSio2-bSiO2*tanh(-cSiO2*ca+dSiO2)
MaxSio2=aSio2-bSiO2*tanh(-cSiO2*1.0+dSiO2)
MinSio2=aSio2-bSiO2*tanh(-cSiO2*0.0+dSiO2)

ndsio=(Sio2-MinSio2)/(MaxSio2-MinSio2)

end function

!!!!!!!=========================================================================
subroutine pdvar(z1,z2,dn,dcn,A1U,B1U,C1U,A1L,B1L,C1L,A1M, B1M, C1M,A1,B1,C1)

implicit none 

integer i,n,z1,z2,dn,dcn, crustn, pd3flag
parameter (n=1000000)
real(8) A1U, B1U, C1U
real(8) A1L, B1L, C1L
real(8) A1M, B1M, C1M

real(8) A1(-n:n)
real(8) B1(-n:n)
real(8) C1(-n:n)
real(8) gradC1

!gradC1=(C1L-C1U)/(abs(crustn-dn))


    
do i=z1,z2
    
    if (i.gt.dn) then
        A1(i)=A1U
        B1(i)=B1U
        C1(i)=C1U
    else if ((i.le.dn).AND.(i.gt.dcn)) then
        A1(i)=A1L
        B1(i)=B1L
        C1(i)=C1L!-gradC1*(i-crustn)
    else
        A1(i)=A1M
        B1(i)=B1M
        C1(i)=C1M
    end if
    
    
end do



        

end subroutine
    
    
    !========================================================
subroutine S(z1,z2,dn,dcn,TsU,TsL,TsM,Ts)
implicit none
real(8) TsU,TsL,TsM
integer i, n,z1,z2,dn,dcn
parameter(n=1000000)
real(8) Ts(-n:n)
real(8)gradTs

!gradTs=(TsL-TsU)/(abs(crustn-dn))


    
do i=z1,z2
    
    if (i.gt.dn) then
        Ts(i)=TsU
    else if ((i.le.dn).AND.(i.gt.dcn)) then
        Ts(i)=TsL!-gradTs*(i-crustn)
    else
        Ts(i)=TsM
    end if 
    
end do


end subroutine
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++

