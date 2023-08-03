import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.constants import kpc
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
mpl.rcParams['axes.linewidth']  = 1
plt.style.use('seaborn-v0_8-deep')
import argparse

parser = argparse.ArgumentParser(description='SED plotting tool')
parser.add_argument('--infile1' , type=str, help='Input (SED) csv file name generated with VOU-Blazars (d/f=--)', default='--')
parser.add_argument('--infile2' , type=str, help='Input (SED) file name formatted for SSDC SED tool (e.g. sed_pow_4SED.txt)', default='--')
parser.add_argument('--infile3' , type=str, help='Input (SED) file from Browse (d/f=None)', default='--')
parser.add_argument('--infile4' , type=str, help='Input (SED) file from UVOT or ASASSN (d/f= None )', default='--')
parser.add_argument('--infile5' , type=str, help='Input average SED file for comparison (d/f= None )', default='--')
parser.add_argument('--infile6' , type=str, help='Input NuSTAR SED (d/f= None )', default='--')
parser.add_argument('--infile7' , type=str, help='Input export SSDC SED tool file (d/f= None )', default='--')
parser.add_argument('--outfile', type=str, help='Sed output file name (d/f=PySED.png)', default='PySED.png')
parser.add_argument('--dataoutfile', type=str, help='Data output file name (d/f=PySED_data.txt)', default='PySED_data.txt')
parser.add_argument('--xaxis'  , type=str, help='X-axis type (f=Frequency,e=Energy eV,k=Energy KeV,g=Energy GeV, t=Energy TeV, d/f=f)', default='f')
parser.add_argument('--title'  , type=str, help='SED title (centred) ', default=' ')
parser.add_argument('--ltitle'  , type=str, help='SED title (left justified)', default=' ')
parser.add_argument('--rtitle'  , type=str, help='SED title (rigth justified)', default=' ')
parser.add_argument('--redshift'  , type=float, help='Redshift', default='0.0')
parser.add_argument('--scaling_factor'  , type=float, help='Scaling factor for data in infil5', default='1.0')
parser.add_argument('--upperlimits'  , type=str, help='Upper limits', default='yes')

parser.add_argument('--GETemplate' , type=str, help='Add a giant elliptical template? (d/f=no)', default='no')
parser.add_argument('--BBTemplate' , type=str, help='Add a Blue Bump + accretion template? (d/f=no)', default='no')
parser.add_argument('--BBScalFactor' , type=float, help='Blue Bump Scaling factor (d/f=5.e-12)', default='5.e-12')
parser.add_argument('--C279Template' , type=str, help='Add a 3C279 template? (LBL, <nupeak> ~1.e13 d/f=no)', default='no')
parser.add_argument('--C279ScalFactor' , type=float, help='3C279 Scaling factor (d/f=1)', default='1')
parser.add_argument('--OJ287Template' , type=str, help='Add a OJ287 template? (IBL <nupeak> ~7.e13 d/f=no)', default='no')
parser.add_argument('--OJ287ScalFactor' , type=float, help='OJ287 Scaling factor (d/f=1)', default='1')
parser.add_argument('--BLLACTemplate' , type=str, help='Add a BL Lac template? (IBL <nupeak> ~1.e14 d/f=no)', default='no')
parser.add_argument('--BLLACScalFactor' , type=float, help='BL Lac Scaling factor (d/f=1)', default='1')
parser.add_argument('--C371Template' , type=str, help='Add a 3C371 template? (IBL <nupeak> ~2.e14 d/f=no)', default='no')
parser.add_argument('--C371ScalFactor' , type=float, help='3C371 Scaling factor (d/f=1)', default='1')
parser.add_argument('--TXS0506Template' , type=str, help='Add a TXS0506+056 template? (IBL, <nupeak> ~5.e14 d/f=no)', default='no')
parser.add_argument('--TXS0506ScalFactor' , type=float, help='TXS0506+056 Scaling factor (d/f=1)', default='1')
parser.add_argument('--PKS2155Template' , type=str, help='Add a PKS2155-304 template? (HBL, <nupeak> ~5.E15 d/f=no)', default='no')
parser.add_argument('--PKS2155ScalFactor' , type=float, help='PKS2155-304 Scaling factor (d/f=1)', default='1')
parser.add_argument('--MKN421Template' , type=str, help='Add a MKN421 template? (HBL, <nupeak> ~1.e17 d/f=no)', default='no')
parser.add_argument('--MKN421ScalFactor' , type=float, help='MKN421 Scaling factor (d/f=1)', default='1')
parser.add_argument('--ES1959Template' , type=str, help='Add a 1ES1959+650 template? (HBL, <nupeak> ~5.e17 d/f=no)', default='no')
parser.add_argument('--ES1959ScalFactor' , type=float, help='1ES1959+650 Scaling factor (d/f=1)', default='1')
parser.add_argument('--MKN501Template' , type=str, help='Add a MRK501 template? (HBL, <nupeak> ~5.e17 d/f=no)', default='no')
parser.add_argument('--MKN501ScalFactor' , type=float, help='MRK501 Scaling factor (d/f=1)', default='1')
parser.add_argument('--ShowLines' , type=str, help='PLot lines at 1e12, 1e13, 1e14 and 1e15 Hz (d/f=no)', default='no')

parser.add_argument('--eblcorrected'  , type=str, help='Show data ebl-corrected', default='no')
parser.add_argument('--erosita'  , type=str, help='Highligh eROSITA data', default='no')
parser.add_argument('--alma'  , type=str, help='Highligh ALMA data', default='no')
parser.add_argument('--smarts'  , type=str, help='Highligh SMARTS data', default='no')
parser.add_argument('--swift'  , type=str, help='Highligh Swift-XRT data', default='no')
parser.add_argument('--xmm'  , type=str, help='Highligh XMM data', default='no')
parser.add_argument('--nustar'  , type=str, help='Highligh NuSTAR data', default='no')

args = parser.parse_args()
type = args.xaxis
infile1 = args.infile1
infile2 = args.infile2
infile3 = args.infile3
infile4 = args.infile4
infile5 = args.infile5
infile6 = args.infile6
infile7 = args.infile7
outfile = args.outfile
dataoutfile = args.dataoutfile
sed_title = args.title
sed_ltitle = args.ltitle
sed_rtitle = args.rtitle
z = args.redshift
scaling_factor = args.scaling_factor
ulimit = args.upperlimits
EBL = args.eblcorrected
eROSITA = args.erosita
ALMA = args.alma
XMM = args.xmm
SWIFT = args.swift
NUSTAR = args.nustar
SMARTS = args.smarts

Templ_GE = args.GETemplate
Templ_BB = args.BBTemplate
ScalingFactorBB = args.BBScalFactor
Templ_3C279 = args.C279Template
ScalingFactor3C279 = args.C279ScalFactor
Templ_3C371 = args.C371Template
ScalingFactor3C371 = args.C371ScalFactor
Templ_OJ287 = args.OJ287Template
ScalingFactorOJ287 = args.OJ287ScalFactor
Templ_BLLAC = args.BLLACTemplate
ScalingFactorBLLAC = args.BLLACScalFactor
Templ_TXS0506 = args.TXS0506Template
ScalingFactorTXS0506 = args.TXS0506ScalFactor
Templ_PKS2155 = args.PKS2155Template
ScalingFactorPKS2155 = args.PKS2155ScalFactor
Templ_MKN421 = args.MKN421Template
ScalingFactorMKN421 = args.MKN421ScalFactor
Templ_MKN501 = args.MKN501Template
ScalingFactorMKN501 = args.MKN501ScalFactor
Templ_1ES1959 = args.ES1959Template
ScalingFactor1ES1959 = args.ES1959ScalFactor
ShowLines = args.ShowLines

f = open(dataoutfile,"w")

if infile1 != '--':
   data  = pd.read_csv(infile1,  delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end", "IsDet", "Cat", "Reference"],  skiprows=1, header=None)
if infile2 != '--':
   data2 = pd.read_csv(infile2,  delimiter="|" , names=["ra","dec","Frequency", "Freq_err", "Flux", "Flux_err", "MJD", "MJD_end","IsDet","a","b"],  skiprows=0, header=None)
if infile3 != '--':
   data3 = pd.read_csv(infile3,  delimiter="," , names=["ra","dec","MJD","F5kev","F5kev_err","F05kev","F05kev_err","F15kev","F15kev_err","F3kev","F3kev_err","F45kev","F45kev_err","F1kev","F1kev_err","IsDet"],  skiprows=0, header=None, index_col=False)
#   data3["Freq_1kev"] = data3["ra"] * 0.0 + 2.418*10**(17)
   data3["Freq_1kev"] =  2.418*10**(17)
   data3["Freq_05kev"] = data3["Freq_1kev"] * 0.5 
   data3["Freq_15kev"] = data3["Freq_1kev"] * 1.5 
   data3["Freq_3kev"] = data3["Freq_1kev"] * 3.0 
   data3["Freq_45kev"] = data3["Freq_1kev"] * 4.5 
if infile4 != '--':
   data4  = pd.read_csv(infile4,  delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end" , "IsDet" , "b"],  skiprows=1, header=None)
if infile5 != '--':
   data5  = pd.read_csv(infile5,  delim_whitespace=True , names=["LogFreq", "zero", "LogFlux", "zero1"],  skiprows=1, header=None)
   data5["xtype"] = 10**data5["LogFreq"]
   data5["Flux"] = 10**data5["LogFlux"]/scaling_factor
if infile6 != '--':
   data6 = pd.read_csv(infile6,  delimiter="|" , names=["ra","dec","Frequency", "Freq_error", "Flux", "Flux_err", "MJD", "d"],  skiprows=0, header=None)
if infile7 != '--':
   data7 = pd.read_csv(infile7, delim_whitespace=True , names=["Frequency", "Ferr", "Flux", "Flux_err", "MJD", "MJD_end"],  skiprows=0, header=None)
fig,ax1 = plt.subplots(figsize = (12, 8), facecolor = 'white', dpi=500)
ax1.set_ylabel(r'E$\cdot$F$_{\rmE}$ [erg cm$^{-2}$ s$^{-1}$]', fontsize=24, fontweight='bold')

if type == 'f':
   if infile1 != '--':
      data["xtype"] = data["Frequency"] 
   if infile2 != '--':
      data2["xtype"] = data2["Frequency"] 
   if infile4 != '--':
      data4["xtype"] = data4["Frequency"] 
   if infile6 != '--':
      data6["xtype"] = data6["Frequency"] 
      dd6 = data6[(data6["Flux"]/data6["Flux_err"] > 2.0) & ( data6["Frequency"] < 1.7e19) ]
   if infile7 != '--':
      data7["xtype"] = data7["Frequency"] 
#   ax1.set_xlabel(r'Frequency ($\nu$, observer frame) [Hz]', fontsize=24, fontweight='bold')
   ax1.set_xlabel(r'$\nu_{\rm\bf observer~frame}{\bf [Hz]}$', fontsize=24 )
   ax1.set_ylabel(r'$\nu$F$_{\nu}$ [erg cm$^{-2}$ s$^{-1}$]', fontsize=24, fontweight='bold')
elif type == 'e':
   data["xtype"] = data["Frequency"] / (2.418*10**(14))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(14))
   ax1.set_xlabel("Energy [eV]" , fontsize=24, fontweight='bold')
elif type == 'k':
   data["xtype"] = data["Frequency"] / (2.418*10**(17))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(17))
   ax1.set_xlabel("Energy [KeV]" , fontsize=24, fontweight='bold')
elif type == 'g':
   data["xtype"] = data["Frequency"] / (2.418*10**(23))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(23))
   ax1.set_xlabel("Energy [GeV]" , fontsize=24, fontweight='bold')
elif type == 't':
   data["xtype"] = data["Frequency"] / (2.418*10**(26))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(26))
   ax1.set_xlabel("Energy [TeV]" , fontsize=24, fontweight='bold')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.tick_params(axis='both', which='both', direction='in', length=5 , labelsize=18, top=True, right=True, labeltop=False, labelright=False)

if infile1 != '--':
   dda = data[data["Cat"] != 'DEBL          ']
   debl = data[data["Cat"] == 'DEBL          ']
   deros = data[(data["Cat"]   == 'eROSITA-EDR   ') & (data["IsDet"]   != 'UL')]
   dalma = data[data["Cat"]   == 'ALMA          ']
   dsmarts = data[data["Cat"] == 'SMARTS        ']
   dnustar = data[data["Cat"] == 'NuBlazar      ']
   dxmm = data[(data["Cat"]  == 'XMMSL2        ') | (data["Cat"]  == '4XMM-DR11     ') ]
   dswift = data[(data["Cat"]  == '1OUSX         ') | (data["Cat"]  == '2SXPS         ')  | (data["Cat"]  == 'OUSXB         ') | (data["Cat"]  == 'XRTSPEC       ')]
   #dd  = dda[(dda["IsDet"] != 'UL') & (dda["Flux_err"] > 0.) ]
   #dd  = dda[( ( (dda["IsDet"] != 'UL') & (dda["Flux_err"] > 0.) ) | ( ( (dda["Frequency"] > 1.36e13) & (dda["Frequency"] < 1.37e13) ) & (dda["Flux"] > 1.e-12) & (dda["IsDet"] != 'UL') ) )]
   dd  = dda[(dda["IsDet"] != 'UL') & (dda["Flux_err"] >= 0.) & ( ( (dda["Frequency"] > 1.4e13) | (dda["Frequency"] < 1.3e13) ) | ( ( (dda["Frequency"] > 1.3e13) & (dda["Frequency"] < 1.4e13) ) & (dda["Flux"] > 1.5e-12) ) ) ]
   dlimit = data[data["IsDet"] == 'UL']
f.write ("Frequency (Hz), nuFnu flux (erg/cm2/s), nuFnu flux error (erg/cm2/s) , MJD (days) \n")
if infile7 != '--':
   ax1.errorbar(data7["xtype"],data7["Flux"],xerr=None,yerr=data7["Flux_err"], fmt='o',color = '#58aaee', markersize='3')
   f.write ("#Data from SSDC SED tool \n")
   df_list = [data7[['xtype', 'Flux', 'Flux_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if infile2 != '--':
   dd2 = data2[(data2["Flux"]/data2["Flux_err"] > 1.8) & ( (data2["Frequency"] < 1.9e18) | (data2["Frequency"] > 1.e21) )]
if infile3 != '--':
   dd3  = data3[(data3["F1kev"]/data3["F1kev_err"] > 1.5) & (data3["IsDet"] != 'UL')]
   d305 = dd3[dd3["F05kev_err"] > 0.]
   d315 = dd3[dd3["F15kev_err"] > 0.]
   d330 = dd3[dd3["F3kev_err"] > 0.]
   d345 = dd3[dd3["F45kev_err"] > 0.]
   dlimit3 = data3[data3["IsDet"] == 'UL']
if infile4 != '--':
   d4  = data4[(data4["IsDet"] != 'UL') & (data4["Flux"]/data4["Flux_err"] > 2.0)]
   dlimit4  = data4[(data4["IsDet"] == 'UL')]
if infile1 != '--':
   ax1.errorbar(dd["xtype"],dd["Flux"],xerr=None,yerr=dd["Flux_err"], fmt='o',color = '#0066cc', markeredgecolor='black', markersize='5')
   f.write ("#Data from VOU-BLazars tool \n")
   df_list = [dd[['xtype', 'Flux', 'Flux_err' , 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   if EBL == 'yes':
      ax1.errorbar(debl["xtype"],debl["Flux"],xerr=None, yerr=debl["Flux_err"], fmt='D',color = '#ff4d4d',fillstyle='none', markersize='5') 
   if eROSITA == 'yes':
      ax1.errorbar(deros["xtype"],deros["Flux"],xerr=None, yerr=deros["Flux_err"], fmt='o',color = '#ff2222', markersize='5') 
   if ALMA == 'yes':
      ax1.errorbar(dalma["xtype"],dalma["Flux"],xerr=None, yerr=dalma["Flux_err"], fmt='o',color = '#009933', markersize='5') 
   if SMARTS == 'yes':
      ax1.errorbar(dsmarts["xtype"],dsmarts["Flux"],xerr=None, yerr=dsmarts["Flux_err"], fmt='o',color = '#00cc99', markersize='5') 
   if SWIFT == 'yes':
      ax1.errorbar(dswift["xtype"],dswift["Flux"],xerr=None, yerr=dswift["Flux_err"], fmt='o',color = '#00ffff', markersize='5') 
   if NUSTAR == 'yes':
      ax1.errorbar(dnustar["xtype"],dnustar["Flux"],xerr=None, yerr=dnustar["Flux_err"], fmt='o',color = '#ff9933', markersize='5') 
   if XMM == 'yes':
      ax1.errorbar(dxmm["xtype"],dxmm["Flux"],xerr=None, yerr=dxmm["Flux_err"], fmt='o',color = '#ff6699', markersize='5') 
if infile2 != '--':
   ax1.errorbar(dd2["xtype"],dd2["Flux"],xerr=None,yerr=dd2["Flux_err"], fmt='o',color = '#1f7604', markersize='3')
   f.write ("#Data from Swift_xrtproc XRT analysis \n")
   df_list = [dd2[['xtype', 'Flux', 'Flux_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f,mode = 'a' , index=None, header=None , sep=',',  float_format='%.8E')
if infile3 != '--':
   ax1.errorbar(dd3["Freq_1kev"] ,dd3["F1kev"] ,xerr=None,yerr=dd3["F1kev_err"] , fmt='o', color = '#fbd799', markersize='4')
   ax1.errorbar(d305["Freq_05kev"],d305["F05kev"],xerr=None,yerr=d305["F05kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   ax1.errorbar(d315["Freq_15kev"],d315["F15kev"],xerr=None,yerr=d315["F15kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   ax1.errorbar(d345["Freq_45kev"],d345["F45kev"],xerr=None,yerr=d345["F45kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   f.write ("#1KeV data from Swift_xrtproc useful for lightcurve use \n")
   df_list = [dd3[['Freq_1kev', 'F1kev', 'F1kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_1kev'] = df_list['Freq_1kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   f.write ("#Data from Swift_xrtproc XIMAGE analysis \n")
   df_list = [d305[['Freq_05kev', 'F05kev', 'F05kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_05kev'] = df_list['Freq_05kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   df_list = [d315[['Freq_15kev', 'F15kev', 'F15kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_15kev'] = df_list['Freq_15kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   df_list = [d345[['Freq_45kev', 'F45kev', 'F45kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_45kev'] = df_list['Freq_45kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None , float_format='%.8E')
ax1.set_xlim(1e7, 5e27)
if infile4 != '--':
   ax1.errorbar(d4["xtype"],d4["Flux"],xerr=None,yerr=d4["Flux_err"], fmt='o',color = '#ffcccc', markersize='5')
#   ax1.errorbar(data4["xtype"],data4["Flux"],xerr=None,yerr=data4["Flux_err"], fmt='o',color = '#ffcccc', markersize='5')
   f.write ("#Data from UVOT or ASASSN  \n")
   df_list = [d4[['xtype', 'Flux', 'Flux_err' , 'MJD']]]
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if infile5 != '--':
   ax1.errorbar(data5["xtype"],data5["Flux"],xerr=None,yerr=None, fmt='o',color = '#ffb84d', markersize='4')
if infile6 != '--':
   ax1.errorbar(dd6["xtype"],dd6["Flux"],xerr=None, yerr=dd6["Flux_err"], fmt='o',color = '#ff00ee', markeredgecolor='black', markersize='7')
   f.write ("#Data from NuSTAR analysis (Middei et al. 2021) \n")
   df_list = [dd6[['xtype', 'Flux', 'Flux_err' , 'MJD']]]
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if ulimit == 'yes':
   if infile1 != '--':
      ax1.errorbar(dlimit["xtype"],dlimit["Flux"], yerr=0.15*dlimit["Flux"], fmt='o', markersize='0', color='#8b908e', uplims=True) 
      f.write ("#Upper limits from VOU_blazars \n")
      f.write ("#Frequency, nufnu flux (erg/cm2/s), MJD (days) \n")
      df_list = [dlimit[['xtype', 'Flux', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['xtype'] = df_list['xtype'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   if infile3 != '--':
      ax1.errorbar(dlimit3["Freq_1kev"],dlimit3["F1kev"],xerr=None,yerr=0.15*dlimit3["F1kev"], fmt='o', markersize='0', color='#888888', uplims=True)
      f.write ("#Upper limits from XRT XIMAGE analysis \n")
      df_list = [dlimit3[['Freq_1kev', 'F1kev', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['Freq_1kev'] = df_list['Freq_1kev'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   if infile4 != '--':
      ax1.errorbar(dlimit4["xtype"],dlimit4["Flux"],xerr=None,yerr=0.15*dlimit4["Flux"], fmt='o', markersize='0', color='#888888', uplims=True)
      f.write ("#Upper limits from ASAS-SN data \n")
      df_list = [dlimit4[['xtype', 'Flux', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['xtype'] = df_list['xtype'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
mnl, mxl = ax1.get_ylim()
f.write ("#Data from template(s) \n")
f.write ("#Frequency, nufnu flux (erg/cm2/s) \n")
if Templ_BB =='yes': 
# alpha_ox = -0.137 Log(nuLnu @1.2e15Hz)+4.704 
#Steffen ApJ 2006, 131,2826
   infile ='/Users/paologiommi/app/Templates/BlueBumpTemplNormalised_4py.txt'
   dataTempl_BB = pd.read_csv(infile, delimiter="\s+", names=["Frequency", "Flux"],  skiprows=1, header=None)
   dataTempl_BB["FreqObserverframe"] = dataTempl_BB["Frequency"] / (1.+z) 
   dataTempl_BB["Flux5000"] = dataTempl_BB["Flux"] * ScalingFactorBB
   if z > 0. :
     distance = cosmo.luminosity_distance(z)*kpc*1e5
     nulnuat1p2e15 = 4.*np.pi*distance.value**2*ScalingFactorBB*1.938
# 1.e15 = 2500 A
# 1.938 is the flux value at nu = 1.2e15 in the BlueBumpTemplNormalised_4py.txt file 
# 0.154 is the flux value at nu = 2.41E17 in the BlueBumpTemplNormalised_4py.txt file
     scalingXflux = 10.**((-0.137*np.log10(nulnuat1p2e15)+4.704+1.0)*2.605)/0.154
     ddUVOTTB = dataTempl_BB [dataTempl_BB["Frequency"] < 1.e16]
     ddXTB = dataTempl_BB [dataTempl_BB["Frequency"] > 1.e16]
     ax1.errorbar(ddUVOTTB["FreqObserverframe"],ddUVOTTB["Flux5000"],xerr=None, yerr=None, fmt='o', color = '#009900', markersize='1') 
     ax1.errorbar(ddXTB["FreqObserverframe"],ddXTB["Flux5000"]*scalingXflux,xerr=None, yerr=None, fmt='o', color = '#990099', markersize='1') 
   ax1.errorbar(dataTempl_BB["FreqObserverframe"],dataTempl_BB["Flux5000"],xerr=None, yerr=None, fmt='o', color = '#009900', markersize='1') 
if Templ_MKN501 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501SyncLow.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
     distance = cosmo.luminosity_distance(z)*kpc*1e5
     distance_mkn501 = cosmo.luminosity_distance(0.033)*kpc*1e5
     scal = -(distance.value/distance_mkn501.value)**-2*ScalingFactorMKN501
     dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"]*scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -1st part\n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501SyncHigh.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal 
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -2nd part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501ICLow.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -3rd part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501ICHigh.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -4th part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_PKS2155 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_OJ287 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_BLLAC =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacSyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacSyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_1ES1959 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template1ES1959SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_MKN421 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_TXS0506 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506SyncLow.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -1st part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506SyncHigh.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -2nd part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506ICLow.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -3rd part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506ICHigh.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -4th part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_3C371 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template3C371SyncLow.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -1st part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371SyncHigh.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -2nd part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371ICLow.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -3rd part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371ICHigh.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -4th part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_3C279 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template3C279SyncLow.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -1st part\n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279SyncHigh.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -2nd part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279ICLow.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -3rd part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279ICHigh.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -4th part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if z != 0.0:
   distance =  cosmo.luminosity_distance(z)*kpc*1e5
   const = 4.*np.pi*distance.value**2
   if Templ_GE =='yes':
      infile ='/Users/paologiommi/app/Templates/GiantEllipticalTemplate_4py.txt'
      dataTempl_GE = pd.read_csv(infile, delimiter="\s+", names=["Frequency", "Luminosity"],  skiprows=1, header=None)
      dataTempl_GE["FreqObserverframe"] = dataTempl_GE["Frequency"] / (1.+z)
      dataTempl_GE["Flux"] = dataTempl_GE["Luminosity"] / const * 1.0
      ax1.errorbar(dataTempl_GE["FreqObserverframe"],dataTempl_GE["Flux"],xerr=None, yerr=None, fmt='o', color = '#ff33cc', markersize='1', ls ='-')
   ax2 = ax1.twinx()
   mn, mx = ax1.get_ylim()
   ax2.set_ylim(mn*const, mx*const)
   ax2.set_yscale('log')
   ax1.tick_params(axis='both', which='both', direction='in', labelsize=18, top=True, right=False, labeltop=False, labelright=False)
   ax2.tick_params(axis='both', which='both', direction='in', labelsize=18, top=True, right=True, labeltop=False )
   if type == 'f':
      ax2.set_ylabel(r'$\nu$L$_{\nu}$ [erg s$^{-1}$]', fontsize=18, fontweight='bold') 
   else:
      ax2.set_ylabel(r'E$\cdot$L$_{\rmE}$ [erg s$^{-1}$]', fontsize=18, fontweight='bold')
mnl, mxl = ax1.get_ylim()

yye = 0.02
if ShowLines =='yes': 
   ax1.axvline(1.e12,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa')
   ax1.axvline(1.e13,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa')
   ax1.axvline(1.e14,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa')
   ax1.axvline(1.e15,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa') 
   ax1.axvline(1.e16,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa') 
   ax1.axvline(1.e17,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa') 
   ax1.axvline(1.e18,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa') 
   yt = mnl*(mxl/mnl)**0.01
   ax1.text(1.e12,yt,'10$^{12}$',fontsize=10 , color = '#888888')
   ax1.text(1.e14,yt,'10$^{14}$',fontsize=10 , color = '#888888')
   ax1.text(1.e15,yt,'10$^{15}$',fontsize=10 , color = '#888888')
   ax1.text(1.e17,yt,'10$^{17}$',fontsize=10 , color = '#888888')
   ax1.text(1.e18,yt,'10$^{18}$',fontsize=10 , color = '#888888')
if (Templ_GE =='yes') & (z != 0.0) : 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   ax1.text(2.e22,yy,'Giant Elliptical template',fontsize=10 , color = '#ff33cc')
if Templ_BB =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   ax1.text(2.e22,yy,'QSO template',fontsize=10 , color = '#009900')
if Templ_MKN501 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='MRK501 template * '+ str(ScalingFactorMKN501)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#ff9900')
if Templ_PKS2155 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='PKS2155-304 template * '+ str(ScalingFactorPKS2155)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#ff0000')
if Templ_TXS0506 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='TXS0506+056 template * '+ str(ScalingFactorTXS0506)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#00e6e6')
if Templ_3C371 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='3C371 template * '+ str(ScalingFactor3C371)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#e6ac00')
if Templ_3C279 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='3C279 template * '+ str(ScalingFactor3C279)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#cc3300')
if Templ_MKN421 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='MKN421 template * '+ str(ScalingFactorMKN421)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#006600')
if Templ_1ES1959 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='1ES1959+650 template * '+ str(ScalingFactor1ES1959)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#0099cc')
if Templ_OJ287 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='OJ287 template * '+ str(ScalingFactorOJ287)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#996633')
if Templ_BLLAC =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='BL Lac template * '+ str(ScalingFactorBLLAC)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#996633')

#plt.title(sed_title, loc='left' , color='#800000', weight='bold')
plt.savefig(outfile, bbox_inches='tight', format='png')

plt.title(sed_title, fontsize=18, fontweight='bold')
plt.title(sed_ltitle, loc='left' , color='#800000', weight='bold')
plt.title(sed_rtitle, loc='right' , color='#800000', weight='bold')
#plt.savefig(outfile, bbox_inches='tight', format='pdf')
plt.savefig(outfile, bbox_inches='tight', format='png')
