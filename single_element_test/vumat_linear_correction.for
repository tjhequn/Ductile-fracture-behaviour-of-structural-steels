C----------------------------------------------------------------------
C     This umat is used for considering the influence of the 
C     stress state on the yield and fracture behaviour of 
C     construction steel!
C     The stress state is characterized in the principal stress and 
C     represented by stress triaxility and lode angle!
C                 written by HE Qun
C                   2020-09-30
c----------------------------------------------------------------------
      subroutine vumat (
c read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     steptime, totaltime, dt, cmname, coordmp, charlength,
     *     props, density, strainInc, relspininc,
     *     tempold, stretchold, defgradold, fieldold,
     *     stressOld, stateOld, enerinternOld, enerinelasOld,
     *     tempnew, stretchnew, defgradnew, fieldnew,
c write only -
     *     stressNew, stateNew, enerinternNew, enerinelasNew )
c
      include 'vaba_param.inc'
c
      dimension coordmp(nblock,*), charlength(nblock), props(nprops),
     1     density(nblock), strainInc(nblock,ndir+nshr),
     2     relspininc(nblock,nshr), tempold(nblock),
     3     stretchold(nblock,ndir+nshr), 
     4     defgradold(nblock,ndir+nshr+nshr),
     5     fieldold(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerinternOld(nblock),
     7     enerinelasOld(nblock), tempnew(nblock),
     8     stretchnew(nblock,ndir+nshr),
     9     defgradnew(nblock,ndir+nshr+nshr),
     1     fieldnew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerinternNew(nblock), enerinelasNew(nblock)
     
C --------------------------------------------------------------------72
      real*8 :: emod, enu, eg, elam, eqplas, eta_n, theta_n, fStrain_n,
     1          omegaD, eleDel, I1_trial, J2_trial, J3_trial, eta_trial, 
     2          theta_trial, rp_trial, drdp_trial, phi_trial, dphideta_trial, 
     3          psi_trial, dpsidtheta_trial, r_trial, fn1_trial, fStrain_n1,
     4          dlambda,fn1, len_ss, eqplas_n1, I1, J2, J3, eta_n1, theta_n1,
     5          rp_n1, drpdp_n1, phi_n1, dphideta_n1, psi_n1, dpsidtheta_n1, 
     6          r_n1, drdlambda, dfn1, pStrain(50), tStress(50), ddsdde(6,6),
     7          matP1(6,6), matP2(6,6), eelas(6), eplas(6), s_trial(6),
     8          ss_trial(6), n_n1(6), s_n1(6), ss_n1(6), deplas(6)
    
      integer :: numEta, numTheta,numData,num_newton, num_test
C     
      character*80 cmname
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,
     1     six=6.0d0,pi=3.1415926d0,newton=15,toler1=1.0d-8,
     2     toler2=1.0d-7,check0 = 1.0d-7)

C----------------------------------------------------------------------
C material constants and parameters definition
C props(1)~props(8):  E,v,and fracture parameters
C props(9)~props(16): number of parameters & triaxility parameters
C props(17)~props(24): number of parameters and lode angle parametes
C props(25)~props(25+2N): number of data point,and 
C                         plastic strain & true stress
C----------------------------------------------------------------------
C material constants and parameters assignment
C---------------------------------------------------------------------- 
      emod = props(1)
      enu  = props(2)
      numEta  = props(9)
      numTheta = props(17)
      numData = props(25)
      data pStrain/50*0.0d0/
      data tStress/50*0.0d0/
      pStrain(1:numData) = props(26:25+numData)
      tStress(1:numData) = props(26+numData:25+2*numData)
C ----------------------------------------------------------------------
C calculate shaer modulus and lamada
C ----------------------------------------------------------------------        
      eg = emod/(two*(one+enu))
      elam = enu*emod/((one+enu)*(one-two*enu))

C----------------------------------------------------------------------
C initialization of essential vector and matrix, such as elstic matrix 
C ddsdde, auxiliary matrix matP1,matP2
C----------------------------------------------------------------------
      do k1=1, 6
          do k2=1, 6
              ddsdde(k2,k1) = zero
              matP1(k2,k1) = zero
              matP2(k2,k1) = zero
          end do
      end do

      do k1=1, 3
          do k2=1, 3
              ddsdde(k2, k1) = elam
              matP2(k2, k1) = -one/three
          end do
          ddsdde(k1, k1) = two*eg+elam
          matP1(k1, k1) = one
          matP2(k1, k1) = two/three
      end do

      do k1=4, 6
          ddsdde(k1, k1) = two*eg
          matP1(k1, k1) = two
          matP2(k1, k1) = one
      end do
C----------------------------------------------------------------------
C initialization of auxiliary constants
C----------------------------------------------------------------------
      root32 = sqrt(three/two)
      root33 = three*sqrt(three)
      root332 = root33/two
C --------------------------------------------------------------------72
C stress update 
C --------------------------------------------------------------------72
      do k = 1, nblock

C----------------------------------------------------------------------
C state variables definition
C statev(1)~statev(6): elstic strain--eelas
C statev(7)~statev(12): plastic strain--eplas
C statev(13): equivalent plastic strain--eqplas
C statev(14): stress triaxility--eta
C statev(15): lode angle--theta
C statev(16): fracture strain--fStrain
C statev(17): fracture indicator--omegaD
C statev(18): elemment deletion--eleDel
C read the results of state variables from last step
C----------------------------------------------------------------------
        eelas = stateOld(k,1:6)
        eplas = stateOld(k,7:12)
        eqplas = stateOld(k,13)
        eta_n = stateOld(k,14)
        theta_n = stateOld(k,15)
        fStrain_n = stateOld(k,16)
        omegaD = stateOld(k,17)
        eleDel = stateOld(k,18)

      if (steptime.eq.zero) then
          eleDel = one
      endif

C calculate trial stress and deviatoric stress

        s_trial = stressOld(k,1:6) + matmul(ddsdde,strainInc(k,1:6))
        ss_trial = matmul(matP2,s_trial)

C calculate stree invariant I1, J2, and J3

        I1_trial = s_trial(1) + s_trial(2) + s_trial(3)
        J2_trial = dot_product(ss_trial,matmul(matP1,ss_trial))/two
        J3_trial = ss_trial(1)*ss_trial(2)*ss_trial(3)+two*
     1             ss_trial(4)*ss_trial(5)*ss_trial(6)-
     2             ss_trial(1)*ss_trial(5)**2-ss_trial(2)*ss_trial(6)**2-
     3             ss_trial(3)*ss_trial(4)**2

C calculate triaxility and lode angle
      
      call stressState(I1_trial,J2_trial,J3_trial,eta_trial,theta_trial)

C calculate the strain hardening 
      call hardening(numData, pStrain, tStress, eqplas, rp_trial, drdp_trial)

C calculate the triaxility effect 
      call phiEta(numEta, props(10), eta_trial, phi_trial,dphideta_trial)

C calculate the lode angle effect
      call psiTheta(numTheta,props(18),theta_trial,psi_trial,dpsidtheta_trial)

C calculate the total hardening effect
      r_trial = rp_trial*phi_trial*psi_trial

C calculate yield
      fn1_trial = sqrt(three*J2_trial) - r_trial
      
C check the trial state
      if (fn1_trial.le.zero) then

C the trial state is the real state, update all state variables

C calculate fracture strian        
        call fractureStrian(props(3), eta_trial, theta_trial, fStrain_n1)

C update all the state variables
        stressNew(k,1:6) = s_trial
        stateNew(k,1:6) = eelas + strainInc(k,1:6)
        stateNew(k,7:12) = eplas
        stateNew(k,13) = eqplas
        stateNew(k,14) = eta_trial
        stateNew(k,15) = theta_trial
        stateNew(k,16) = fStrain_n1
        stateNew(k,17) = omegaD
        stateNew(k,18) = eleDel
        stateNew(k,19) = phi_trial
        stateNew(k,20) = psi_trial

        stressPower = dot_product((stressOld(k,1:6)+
     1  stressNew(k,1:6))/two,strainInc(k,1:6))
        enerinternNew(k) = enerinternOld(k)+ stressPower/density(k)
        enerinelasNew(k) = enerinelasOld(k)

      else
C the trial state is not the real state, perform plastic correction
C initialize essential variables
        num_newton = 0
        dlambda = 0.0d0
        fn1 = 100.0d0

C calculate the direciton of plastic flow
        len_ss = sqrt(two*J2_trial)
        if (len_ss.le.check0) then
          n_n1 = ss_trial*len_ss
        else
          n_n1 = ss_trial/len_ss
        endif

C----------------------------------------------------------------------
C plastic correction
C----------------------------------------------------------------------
        do while ((abs(fn1).gt.toler1).and.(num_newton.le.newton))
          num_newton = num_newton + 1
          if (num_newton.gt.newton) then
            print *,'***error - too many attempts for vumat'
            stop
          endif

C calculate deviatoric stress
          s_n1 = s_trial - sqrt(6.0d0)*eg*dlambda*n_n1
          ss_n1 = matmul(matP2,s_n1)
          eqplas_n1 = eqplas + dlambda

          I1 = s_n1(1) + s_n1(2) + s_n1(3)
          J2 = dot_product(ss_n1,matmul(matP1,ss_n1))/two
          J3 = ss_n1(1)*ss_n1(2)*ss_n1(3)+two*
     1         ss_n1(4)*ss_n1(5)*ss_n1(6)-
     2         ss_n1(1)*ss_n1(5)**2-ss_n1(2)*ss_n1(6)**2-
     3         ss_n1(3)*ss_n1(4)**2

C calculate triaxility and lode angle
      
          call stressState(I1,J2,J3,eta_n1,theta_n1)
          theta_n1 = theta_trial

C calculate the strain hardening
          call hardening(numData,pStrain,tStress,eqplas_n1,rp_n1,drpdp_n1)

C calculate the triaxility effect 
          call phiEta(numEta,props(10),eta_n1,phi_n1,dphideta_n1)

C calculate the lode angle effect
          call psiTheta(numTheta,props(18),theta_n1,psi_n1,dpsidtheta_n1)

C calculate the total hardening effect
          r_n1 = rp_n1*phi_n1*psi_n1

C calculate yield
          fn1 = sqrt(three*J2) - r_n1

C check the update fn1

          if (abs(fn1).gt.toler1) then

C update delta_lambda if fn1 does not stastify the condition
C calculate derivates dfn1

            drdlambda = drpdp_n1*phi_n1*psi_n1+eg*I1_trial*rp_n1*dphideta_n1*psi_n1/3.0d0/J2
            dfn1 = -3.0d0*eg - drdlambda
            dlambda = dlambda - fn1/dfn1

          else

            deplas= root32*dlambda*n_n1
            eplas = eplas + deplas
            eelas = eelas + strainInc(k,1:6) - deplas

            call fractureStrian(props(3), eta_n1, theta_n1, fStrain_n1)
C             write(*,*) 'fStrain_n1 = ', fStrain_n1
            omegaD = omegaD + dlambda*(1.0d0/fStrain_n+1.0d0/fStrain_n1)/2.0d0

            if (omegaD.ge.one) then
              eleDel = zero
            else
              eleDel = one
            endif

            stressNew(k,1:6) = s_n1
            stateNew(k,1:6) = eelas
            stateNew(k,7:12) = eplas
            stateNew(k,13) = eqplas_n1
            stateNew(k,14) = eta_n1
            stateNew(k,15) = theta_n1
            stateNew(k,16) = fStrain_n1
            stateNew(k,17) = omegaD
            stateNew(k,18) = eleDel
            stateNew(k,19) = phi_n1
            stateNew(k,20) = psi_n1

            stressPower = dot_product((stressOld(k,1:6)+
     1                 stressNew(k,1:6)),strainInc(k,1:6))/two
            enerinternNew(k)=enerinternOld(k)+stressPower/density(k)                           
            plasticworkinc = dot_product((stressOld(k,1:6)+
     1                     stressNew(k,1:6)),deplas)/two
            enerinelasNew(k) = enerinelasOld(k)+
     1                            plasticworkinc/density(k)

            endif
          end do
        endif
      end do

      return
      end


C----------------------------------------------------------------------
C  stressState is used to calculate the triaxility and lode angle and 
C  related variables
C----------------------------------------------------------------------
      subroutine stressState(I1,J2,J3,eta,theta)

        include 'vaba_param.inc'

        real*8 I1,J2,J3,eta,theta,cos3theta,root332

        parameter(one=1.0d0,two=2.0d0,three=3.0d0,pi=3.1415926d0,check0=1.0e-10)

        root332 = three*sqrt(three)/two

        if (J2.le.check0) then
          eta = 0.0d0
          theta = 0.0d0
        else
          eta = I1/sqrt(three*J2)/three

          cos3theta = root332*J3/J2**(1.5d0)
          if (abs(cos3theta).gt.1.0d0) then
            cos3theta = sign(1.0d0,cos3theta)
            theta = acos(cos3theta)/three
          else
            theta = acos(cos3theta)/three
          endif
          theta = one - 6.0d0*theta/pi

        endif 

        return
      end
C--------------------------end of stressState--------------------------

C----------------------------------------------------------------------
C  hardening is a subroutine used for calculating the strain hardening 
C  of a given state without considering the effect of stress state
C----------------------------------------------------------------------
      subroutine hardening(N,pstrain,stress,peeq,ryield,drdp)
        
        include 'vaba_param.inc'
        integer N,loc
        real*8 :: pstrain(N),stress(N),peeq,ryield,drdp

        if (peeq.ge.pstrain(N)) then
            ryield = stress(N)
            drdp = 0.d0
        else
            loc  = maxloc(pstrain,dim=1,mask=pstrain.le.peeq)
            drdp = (stress(loc+1)-stress(loc))/
     1               (pstrain(loc+1)-pstrain(loc))
            ryield = stress(loc)+drdp*(peeq-pstrain(loc))
        endif

        return
      end
C---------------------------end of hardening---------------------------



C----------------------------------------------------------------------
C  phiEta is a subroutine used for considering the effect of triaxility 
C----------------------------------------------------------------------
      subroutine phiEta(N,etaPara,eta,phi,dphideta)
        
        include 'vaba_param.inc'

        integer N
        dimension :: etaPara(N)
        real*8 :: eta, eta0, c_eta, phi,dphideta       
C----------------------------------------------------------------------
C Bai model is adopted here for considering the stress triaxility 
C effect. Therefore only two parameters is needed. if you have a new  
C model to examine the effect of triaxility, only the following code  
C needs to be modified for your specific proposes.
C It should be noted here you can only modified this part if the number 
C of your parameters are less than seven.
C----------------------------------------------------------------------
C considering the effect of triaxility
        c_eta = etaPara(1)
        eta0  = etaPara(2)

        phi = 1.0d0 - c_eta*(eta - eta0)
        dphideta = -c_eta

        return
      end
C----------------------------end of phiEta----------------------------- 


C----------------------------------------------------------------------
C  psiTheta is a subroutine used for considering lode angle effect 
C----------------------------------------------------------------------
      subroutine psiTheta(N,thetaPara,theta,psi,dpsidtheta)
        
        include 'vaba_param.inc'

        integer N
        dimension :: thetaPara(N)
        real*8 :: theta, psi, dpsidtheta, gamma,ctheta_s,
     1            ctheta_t, ctheta_c, m, psi_norm, dgammadtheata
        parameter (pi=3.1415926d0, psi_const=6.464101615d0)
C----------------------------------------------------------------------
C Bai model is adopted here for considering the lode angle effect.
C Therefore only four parameters is needed. if you have a new  
C model to examine the effect of lode angle, only the following code  
C needs to be modified for your specific proposes.
C It should be noted here you can only modified this part if the number 
C of your parameters are less than seven.
C----------------------------------------------------------------------
        ctheta_s = thetaPara(1)
        ctheta_t = thetaPara(2)
        ctheta_c = thetaPara(3)
        m = thetaPara(4)

        if (theta.ge.0.0d0) then
          ctheta_ax = ctheta_t
          dpsidtheta = (ctheta_ax-ctheta_s)
        else
          ctheta_ax = ctheta_c
          dpsidtheta = -(ctheta_ax-ctheta_s)
        endif

        psi = ctheta_s+(ctheta_ax-ctheta_s)*abs(theta)

        return
      end
C--------------------------end of psiTheta--------------------------- 


C----------------------------------------------------------------------
C  fractureStrian is a subroutine used for calculate fracture strain 
C----------------------------------------------------------------------
      subroutine fractureStrian(fracturePara,eta,theta,fstrain)
        
        include 'vaba_param.inc'

        dimension :: fracturePara(6)
        real*8 :: A,n,c1,c2
        real*8 :: eta, theta, fstrain
        parameter (pi=3.1415926d0)
C----------------------------------------------------------------------
C Modified Mohr-Coulomb rule is adopted here.Therefore only four   
C parameters is needed. if you have a new model to predict fracture   
C strain, only the following code needs to be modified for your 
C specific proposes. It should be noted here you can only modified this 
C part if the number of your parameters are less than six.
C----------------------------------------------------------------------
        A  = fracturePara(1)
        n  = fracturePara(2)
        c1 = fracturePara(3)
        c2 = fracturePara(4)        
        fstrain=(A*(sqrt((1.0d0+c1**2)/3.0d0)*cos(pi*theta/6.0d0)+
     1           c1*(eta+sin(pi*theta/6.0d0)/3.0d0))/c2)**(-1.0d0/n)
        return
      end
C--------------------------end of fractureStrian-----------------------