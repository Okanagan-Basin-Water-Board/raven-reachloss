"""Apply the OBWB reach-loss edits (D181) to a clean Raven v4.12 source tree.
Usage: python3 apply_reachloss.py <raven_root>   (edits <raven_root>/src in place; run on a COPY)
Every edit must match exactly once, or the script stops (FLAG)."""
import sys, os
ROOT = sys.argv[1]
M = "// OBWB patch D181 (reach loss)"

def sub(rel, old, new):
    path = os.path.join(ROOT, "src", rel)
    t = open(path, newline="").read()
    crlf = "\r\n" in t
    if crlf: old = old.replace("\n", "\r\n"); new = new.replace("\n", "\r\n")
    n = t.count(old)
    if n != 1:
        sys.exit(f"[FLAG] {rel}: expected 1 match, found {n}: {old.strip()[:70]!r}")
    open(path, "w", newline="").write(t.replace(old, new.replace("@M@", M)))
    print(f"[PASS] edit {rel}: {old.strip()[:60]!r}")

# ---------------- SubBasin.h ----------------
sub("SubBasin.h",
"""  diversion()  {aQsource=NULL; aQdivert=NULL;julian_start=julian_end=0; target_p=DOESNT_EXIST;min_flow=0;percentage=1.0;nPoints=0;}""",
"""  @M@: fields below are used only when loss_form>0 (:ReachLoss); legacy diversions keep loss_form=0
  int       loss_form;      //0 legacy diversion; 1 EMPIRICAL_Q; 2 DISCONNECTED; 3 CONNECTED
  double    lp[12];         //loss parameters, see :ReachLoss
  int       dest_type;      //0 SINK; 1 HRUGROUP storage; 2 SUBBASIN via linear return store
  int       dest_group;     //global HRU group index (dest_type 1)
  int       dest_iSV;       //state variable index (dest_type 1)
  int       dest_p;         //target subbasin index (dest_type 2)
  double    dest_k;         //return store time constant [d] (dest_type 2)
  double    store;          //return store volume [m3] (dest_type 2)
  double    Qact,QactLast;  //actual loss rate, this and last step [m3/s]
  double    Qpre;           //start-of-step pre-diversion outflow used by the loss law this step [m3/s]
  double    vol_step;       //volume lost this step [m3], same trapezoid as the channel balance
  double    vol_pending;    //volume awaiting deposit to the HRU group [m3]
  long long nGainClip;      //steps where the CONNECTED head difference was negative (gain clipped to 0)

  diversion()  {aQsource=NULL; aQdivert=NULL;julian_start=julian_end=0; target_p=DOESNT_EXIST;min_flow=0;percentage=1.0;nPoints=0;
                loss_form=0; for(int j=0;j<12;j++){lp[j]=0.0;} dest_type=0; dest_group=DOESNT_EXIST; dest_iSV=DOESNT_EXIST; dest_p=DOESNT_EXIST;
                dest_k=1.0; store=0.0; Qact=QactLast=0.0; Qpre=0.0; vol_step=0.0; vol_pending=0.0; nGainClip=0;}""")
sub("SubBasin.h",
"""///////////////////////////////////////////////////////////////////
/// \\brief flow diversion data strucure""",
"""///////////////////////////////////////////////////////////////////
/// \\brief flow diversion data strucure""")  # anchor check only
sub("SubBasin.h",
"""  double          GetDiversionFlow         (const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &pDivert) const;""",
"""  double          GetDiversionFlow         (const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &pDivert) const;
  @M@
  void            AddReachLoss             (const int jul_start, const int jul_end, const int form, const double *par, const int npar,
                                            const int dest_type, const int dest_group, const int dest_iSV, const int dest_p, const double dest_k);
  diversion      *GetDiversionPtr          (const int i) const;
  bool            HasReachLoss             () const;
  double          GetInternalReachLossVolume() const;   //[m3] lost this step to HRUGROUP or SUBBASIN destinations
  double          GetReachLossStorage      () const;   //[m3] in return stores and pending HRU deposits""")

# ---------------- SubBasin.cpp ----------------
sub("SubBasin.cpp",
"""double CSubBasin::GetDiversionFlow(const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &p_Divert) const
{""",
"""@M@
// Loss rate [m3/s] from the reach outflow estimate Q [m3/s]. Never negative, never above Q (supply limit).
//  EMPIRICAL_Q : lp = a, b, cap                                              L = min(a Q^b, cap)
//  DISCONNECTED: lp = Kbed[m/d], d[m], length[m], aW, bW, aH, bH, cap        W = aW Q^bW, h = aH Q^bH
//                L = Kbed (h + d)/d W length                                 (water table below the bed bottom)
//  CONNECTED   : as DISCONNECTED plus lp[8]=z_bed[m], lp[9]=h_gw[m]
//                L = Kbed ((z_bed + h) - max(h_gw, z_bed - d))/d W length    (continuous with DISCONNECTED)
//  cap <= 0 means no cap. A negative head difference (gaining) is clipped to 0 and flagged through 'gain'.
static double CalcReachLossRate(const diversion *d, const double &Q, bool &gain)
{
  gain=false;
  if (!(Q>0.0)) { return 0.0; }   //zero, negative or NaN flow: no loss
  double Qs=Q;
  double L=0.0;
  if (d->loss_form==1) {
    L=d->lp[0]*pow(Qs,d->lp[1]);
    if (d->lp[2]>0.0) { L=min(L,d->lp[2]); }
  }
  else {
    double K=d->lp[0], b=d->lp[1], len=d->lp[2];
    double W=d->lp[3]*pow(Qs,d->lp[4]);
    double h=d->lp[5]*pow(Qs,d->lp[6]);
    double dh;
    if (d->loss_form==2) { dh=h+b; }
    else                 { dh=(d->lp[8]+h)-max(d->lp[9],d->lp[8]-b); }
    L=K*dh/b*W*len/SEC_PER_DAY;
    if (d->lp[7]>0.0) { L=min(L,d->lp[7]); }
  }
  if (!(L==L)) { ExitGracefully("ReachLoss: loss rate is not a number (check :Parameters)",RUNTIME_ERR); }
  if (L<0.0) { gain=true; L=0.0; }
  return min(L,Qs);
}
double CSubBasin::GetDiversionFlow(const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &p_Divert) const
{""")
sub("SubBasin.cpp",
"""  if (_pDiversions[i]==NULL) { p_Divert=DOESNT_EXIST; return 0.0; } //no diversion

  p_Divert=_pDiversions[i]->target_p;
""",
"""  if (_pDiversions[i]==NULL) { p_Divert=DOESNT_EXIST; return 0.0; } //no diversion

  @M@: reach losses never route directly to a subbasin here; destinations are handled in Solvers.cpp
  if (_pDiversions[i]->loss_form>0) {
    p_Divert=DOESNT_EXIST;
    if (!IsInDateRange(tt.julian_day,_pDiversions[i]->julian_start,_pDiversions[i]->julian_end)) { return 0.0; }
    bool gain;
    return CalcReachLossRate(_pDiversions[i],Q,gain);
  }

  p_Divert=_pDiversions[i]->target_p;
""")
sub("SubBasin.cpp",
"""void  CSubBasin::AddFlowDiversion(const int jul_start,const int jul_end,const int target_p,const double min_flow,const double pct) {
""",
"""@M@
void  CSubBasin::AddReachLoss(const int jul_start, const int jul_end, const int form, const double *par, const int npar,
                              const int dest_type, const int dest_group, const int dest_iSV, const int dest_p, const double dest_k)
{
  ExitGracefullyIf((jul_start<0) || (jul_start>366) || (jul_end<0) || (jul_end>366),"CSubBasin::AddReachLoss: julian dates must be between 0 and 366",BAD_DATA);
  ExitGracefullyIf((form<1) || (form>3),"CSubBasin::AddReachLoss: unknown loss form",BAD_DATA);
  ExitGracefullyIf(npar>12,"CSubBasin::AddReachLoss: too many parameters",BAD_DATA);
  ExitGracefullyIf((dest_type==2) && (dest_k<=0.0),"CSubBasin::AddReachLoss: return store time constant must be positive",BAD_DATA);
  diversion *div=new diversion();
  div->julian_start=jul_start;
  div->julian_end  =jul_end;
  div->target_p    =DOESNT_EXIST;
  div->loss_form   =form;
  for (int j=0;j<npar;j++){ div->lp[j]=par[j]; }
  div->dest_type   =dest_type;
  div->dest_group  =dest_group;
  div->dest_iSV    =dest_iSV;
  div->dest_p      =dest_p;
  div->dest_k      =dest_k;
  if(!DynArrayAppend((void**&)(_pDiversions),(void*)(div),_nDiversions)) {
    ExitGracefully("CSubBasin::AddReachLoss: adding NULL reach loss",BAD_DATA);
  }
}
diversion *CSubBasin::GetDiversionPtr(const int i) const {
  if ((i<0) || (i>=_nDiversions)){ExitGracefully("CSubBasin::GetDiversionPtr: bad index",RUNTIME_ERR);}
  return _pDiversions[i];
}
bool CSubBasin::HasReachLoss() const {
  for (int i=0;i<_nDiversions;i++){ if ((_pDiversions[i]!=NULL) && (_pDiversions[i]->loss_form>0)){return true;} }
  return false;
}
double CSubBasin::GetInternalReachLossVolume() const {
  double sum=0.0;
  for (int i=0;i<_nDiversions;i++){
    if ((_pDiversions[i]!=NULL) && (_pDiversions[i]->loss_form>0) && (_pDiversions[i]->dest_type>0)){ sum+=_pDiversions[i]->vol_step; }
  }
  return sum;
}
double CSubBasin::GetReachLossStorage() const {
  double sum=0.0;
  for (int i=0;i<_nDiversions;i++){
    if ((_pDiversions[i]!=NULL) && (_pDiversions[i]->loss_form>0)){ sum+=_pDiversions[i]->store+_pDiversions[i]->vol_pending; }
  }
  return sum;
}
void  CSubBasin::AddFlowDiversion(const int jul_start,const int jul_end,const int target_p,const double min_flow,const double pct) {
""")
sub("SubBasin.cpp",
"""  //Update flows
  //------------------------------------------------------
  _QoutLast=_aQout[_nSegments-1];""",
"""  @M@: pre-diversion outflow at the start of the step, the Q that Solvers.cpp passed to GetDiversionFlow
  double Qpre_RL=_aQout[_nSegments-1]+_Qdiverted;

  //Update flows
  //------------------------------------------------------
  _QoutLast=_aQout[_nSegments-1];""")
sub("SubBasin.cpp",
"""  _QdivLast=_Qdiverted;
  _Qdiverted=divQ_act;
""",
"""  _QdivLast=_Qdiverted;
  _Qdiverted=divQ_act;

  @M@: book each reach loss's share of the actual (supply-limited) diversion, same trapezoid as the channel balance
  if (HasReachLoss())
  {
    double req_total=0.0;
    int    pd;
    for (int i=0;i<_nDiversions;i++){ req_total+=GetDiversionFlow(i,Qpre_RL,Options,tt,pd); }
    double scale=(req_total>0.0) ? (divQ_act/req_total) : 0.0;
    for (int i=0;i<_nDiversions;i++){
      diversion *d=_pDiversions[i];
      if ((d==NULL) || (d->loss_form==0)){continue;}
      double req=GetDiversionFlow(i,Qpre_RL,Options,tt,pd);
      bool gain=false;
      if ((d->loss_form==3) && IsInDateRange(tt.julian_day,d->julian_start,d->julian_end)){ CalcReachLossRate(d,Qpre_RL,gain); }
      if (gain && !initialize){ d->nGainClip++; }
      d->QactLast=d->Qact;
      d->Qact    =req*scale;
      d->Qpre    =Qpre_RL;
      if (initialize){ d->QactLast=d->Qact; d->vol_step=0.0; continue; }
      d->vol_step=0.5*(d->Qact+d->QactLast)*tstep*SEC_PER_DAY;
      if      (d->dest_type==1){ d->vol_pending+=d->vol_step; }
      else if (d->dest_type==2){ d->store      +=d->vol_step; }
    }
  }
""")

# ---------------- Model.cpp ----------------
sub("Model.cpp",
"""          _CumulOutput-=Qdiv/area*MM_PER_METER; //water that is not diverted out of the watershed
        }
      }""",
"""          _CumulOutput-=Qdiv/area*MM_PER_METER; //water that is not diverted out of the watershed
        }
      }
      @M@: reach losses sent to HRU storage or a return store stay in the watershed
      _CumulOutput-=_pSubBasins[p]->GetInternalReachLossVolume()/area*MM_PER_METER;""")
sub("Model.cpp",
"""    sum+=_pSubBasins[p]->GetChannelStorage();            //[m3]
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;""",
"""    sum+=_pSubBasins[p]->GetChannelStorage();            //[m3]
    sum+=_pSubBasins[p]->GetReachLossStorage();          //[m3] @M@: return stores (0 unless :ReachLoss used)
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;""")

# ---------------- Solvers.cpp ----------------
sub("Solvers.cpp",
"""    //User-specified inflows to upstream end of subbasin reach
    aQinnew[p]+=pBasin->GetSpecifiedInflow(t+tstep);""",
"""    @M@: release from linear return stores (dest SUBBASIN) as lateral inflow volume to the target this step
    for(int i=0; i<pBasin->GetNumDiversions();i++) {
      diversion *dRL=pBasin->GetDiversionPtr(i);
      if ((dRL!=NULL) && (dRL->loss_form>0) && (dRL->dest_type==2) && (dRL->store>0.0)) {
        double R=dRL->store*(1.0-exp(-tstep/dRL->dest_k)); //[m3]
        dRL->store-=R;
        aRouted[dRL->dest_p]+=R;
      }
    }
    //User-specified inflows to upstream end of subbasin reach
    aQinnew[p]+=pBasin->GetSpecifiedInflow(t+tstep);""")
sub("Solvers.cpp",
"""  pModel->AssimilationBackPropagate(Options,tt); //modifies flows using upstream propagation, if needed
""",
"""  pModel->AssimilationBackPropagate(Options,tt); //modifies flows using upstream propagation, if needed

  @M@: deposit losses into HRU-group storage (area weighted over enabled HRUs) and write ReachLoss.csv
  {
    static ofstream RL;
    static bool RLopen=false;
    for (p=0;p<NB;p++) {
      pBasin=pModel->GetSubBasin(p);
      if (!pBasin->HasReachLoss()){continue;}
      for(int i=0; i<pBasin->GetNumDiversions();i++) {
        diversion *dRL=pBasin->GetDiversionPtr(i);
        if ((dRL==NULL) || (dRL->loss_form==0)){continue;}
        if ((dRL->dest_type==1) && (dRL->vol_pending>0.0)) {
          CHRUGroup *pGrp=pModel->GetHRUGroup(dRL->dest_group);
          double Agrp=0.0;
          for (int kk=0;kk<pGrp->GetNumHRUs();kk++){ if (pGrp->GetHRU(kk)->IsEnabled()){ Agrp+=pGrp->GetHRU(kk)->GetArea(); } }
          ExitGracefullyIf(Agrp<=0.0,"ReachLoss: destination HRU group has no enabled area",RUNTIME_ERR);
          for (int kk=0;kk<pGrp->GetNumHRUs();kk++){
            CHydroUnit *pH=pGrp->GetHRU(kk);
            if (pH->IsEnabled()){ aPhinew[pH->GetGlobalIndex()][dRL->dest_iSV]+=dRL->vol_pending/(Agrp*M2_PER_KM2)*MM_PER_METER; }
          }
          dRL->vol_pending=0.0;
        }
        if (!RLopen) {
          if (nConstituents>0) {
            WriteWarning("ReachLoss: water removed by :ReachLoss carries no constituent mass or heat; constituent mass balances will not close",Options.noisy);
          }
          RL.open(FilenamePrepare("ReachLoss.csv",Options).c_str());
          ExitGracefullyIf(RL.fail(),"ReachLoss: cannot open ReachLoss.csv for writing",FILE_OPEN_ERR);
          RL<<"date,SBID,index,form,dest_type,Qpre_m3s,Qact_m3s,vol_step_m3,store_m3,gain_clip_steps"<<endl;
          RLopen=true;
        }
        RL<<tt.date_string<<","<<pBasin->GetID()<<","<<i<<","<<dRL->loss_form<<","<<dRL->dest_type<<","
          <<setprecision(12)<<dRL->Qpre<<","<<dRL->Qact<<","<<dRL->vol_step<<","<<dRL->store<<","<<dRL->nGainClip<<endl;
      }
    }
  }
""")

# ---------------- ParseTimeSeriesFile.cpp ----------------
sub("ParseTimeSeriesFile.cpp",
"""    else if  (!strcmp(s[0],":UserTimeSeries"              )){code=74; }""",
"""    else if  (!strcmp(s[0],":UserTimeSeries"              )){code=74; }
    else if  (!strcmp(s[0],":ReachLoss"                   )){code=75; } // @M@""")
sub("ParseTimeSeriesFile.cpp",
"""    case (69): //---------------------------------------------
    {/*:UnusableFlowPercentage [SBID] [fraction]*/""",
"""    case (75): //--------------------------------------------- @M@
    {/*:ReachLoss [SBID] [EMPIRICAL_Q|DISCONNECTED|CONNECTED] {start_day} {end_day}
        :Parameters [p1] [p2] ...      (3, 8 or 10 values by form; see SubBasin.cpp CalcReachLossRate)
        :Destination SINK
        :Destination HRUGROUP [group name] [state variable, e.g. SOIL[2]]
        :Destination SUBBASIN [SBID] [k, days]
      :EndReachLoss */
      if(Options.noisy) { cout <<"Reach loss"<<endl; }
      ExitGracefullyIf(Len<3,":ReachLoss: needs [SBID] [FORM]",BAD_DATA);
      CSubBasin *pSB=pModel->GetSubBasinByID(s_to_ll(s[1]));
      ExitGracefullyIf(pSB==NULL,":ReachLoss: subbasin not in model",BAD_DATA);
      int form=0;
      if      (!strcmp(s[2],"EMPIRICAL_Q" )){form=1;}
      else if (!strcmp(s[2],"DISCONNECTED")){form=2;}
      else if (!strcmp(s[2],"CONNECTED"   )){form=3;}
      ExitGracefullyIf(form==0,":ReachLoss: form must be EMPIRICAL_Q, DISCONNECTED or CONNECTED",BAD_DATA);
      int rl_start=0, rl_end=366;
      if (Len>=5){ rl_start=s_to_i(s[3]); rl_end=s_to_i(s[4]); }
      double par[12]; int npar=0;
      for (int j=0;j<12;j++){ par[j]=0.0; }
      int dtype=-1, dgrp=DOESNT_EXIST, dsv=DOESNT_EXIST, dp=DOESNT_EXIST; double dk=1.0;
      while (true) {
        bool eof=p->Tokenize(s,Len);
        ExitGracefullyIf(eof,":ReachLoss: missing :EndReachLoss",BAD_DATA);
        if ((Len==0) || IsComment(s[0],Len)){continue;}
        if (!strcmp(s[0],":EndReachLoss")){break;}
        else if (!strcmp(s[0],":Parameters")){
          npar=Len-1;
          ExitGracefullyIf(npar>12,":ReachLoss: too many parameters",BAD_DATA);
          for (int j=0;j<npar;j++){ par[j]=s_to_d(s[j+1]); }
        }
        else if (!strcmp(s[0],":Destination")){
          ExitGracefullyIf(Len<2,":ReachLoss: :Destination needs a type",BAD_DATA);
          if (!strcmp(s[1],"SINK")){ dtype=0; }
          else if (!strcmp(s[1],"HRUGROUP")){
            ExitGracefullyIf(Len<4,":ReachLoss: :Destination HRUGROUP [name] [state variable]",BAD_DATA);
            CHRUGroup *pGrp=pModel->GetHRUGroup(s[2]);
            ExitGracefullyIf(pGrp==NULL,":ReachLoss: destination HRU group not found",BAD_DATA);
            int lay=0;
            sv_type typ=pModel->GetStateVarInfo()->StringToSVType(s[3],lay,true);
            dsv=pModel->GetStateVarIndex(typ,lay);
            ExitGracefullyIf(dsv==DOESNT_EXIST,":ReachLoss: destination state variable not in model",BAD_DATA);
            ExitGracefullyIf(!CStateVariable::IsWaterStorage(typ),":ReachLoss: destination state variable is not a water store",BAD_DATA);
            dtype=1; dgrp=pGrp->GetGlobalIndex();
          }
          else if (!strcmp(s[1],"SUBBASIN")){
            ExitGracefullyIf(Len<4,":ReachLoss: :Destination SUBBASIN [SBID] [k]",BAD_DATA);
            dp=pModel->GetSubBasinIndex(s_to_ll(s[2]));
            ExitGracefullyIf(dp==DOESNT_EXIST,":ReachLoss: destination subbasin not in model",BAD_DATA);
            dk=s_to_d(s[3]);
            ExitGracefullyIf(dk<=0.0,":ReachLoss: return store k must be positive",BAD_DATA);
            dtype=2;
          }
          else { ExitGracefully(":ReachLoss: destination must be SINK, HRUGROUP or SUBBASIN",BAD_DATA); }
        }
        else { ExitGracefully(":ReachLoss: unrecognized line inside block",BAD_DATA); }
      }
      ExitGracefullyIf((rl_start<0) || (rl_start>366) || (rl_end<0) || (rl_end>366),":ReachLoss: julian start and end days must be between 0 and 366",BAD_DATA);
      if (Options.management_optimization) {
        WriteWarning(":ReachLoss has not been tested with management optimization; the loss is applied outside the optimization problem",Options.noisy);
      }
      int need=(form==1) ? 3 : ((form==2) ? 8 : 10);
      ExitGracefullyIf(npar!=need,":ReachLoss: wrong number of :Parameters for this form (EMPIRICAL_Q 3, DISCONNECTED 8, CONNECTED 10)",BAD_DATA);
      ExitGracefullyIf(dtype<0,":ReachLoss: :Destination missing",BAD_DATA);
      for (int j=0;j<npar;j++){ ExitGracefullyIf(!(par[j]==par[j]),":ReachLoss: a parameter is not a number",BAD_DATA); }
      if (form>1){ ExitGracefullyIf((par[0]<0.0) || (par[1]<=0.0) || (par[2]<=0.0) || (par[3]<0.0) || (par[5]<0.0),
                                    ":ReachLoss: Kbed, aW, aH must be >= 0; bed thickness and length must be > 0",BAD_DATA);
                   ExitGracefullyIf((par[4]<0.0) || (par[4]>1.0) || (par[6]<0.0) || (par[6]>1.0),
                                    ":ReachLoss: hydraulic-geometry exponents bW and bH must be between 0 and 1",BAD_DATA); }
      else       { ExitGracefullyIf((par[0]<0.0) || (par[1]<0.0),":ReachLoss: coefficient a and exponent b must be >= 0",BAD_DATA); }
      pSB->AddReachLoss(rl_start,rl_end,form,par,npar,dtype,dgrp,dsv,dp,dk);
      break;
    }
    case (69): //---------------------------------------------
    {/*:UnusableFlowPercentage [SBID] [fraction]*/""")
print("[PASS] all edits applied")
