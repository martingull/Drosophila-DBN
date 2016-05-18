% LoadFlyData.m creates several global variables 
% for use in the computer lab:
% 
% GeneExpr is a 7 by 58 by 7 matrix giving the (normalized) expression
% of the 7 genes (Hb, Kr, Gt, Kni, Bcd, Cad, Tll) at 58 spatial points
% along the A-P axis of the Drosophila embryo (from 35% to 92% of embryo
% length, in 1% increments), and at seven times during early development
% (spanning about an hour, but we won't worry about the true time scale).
%
% GeneNames is a cell array of strings holding the abbreviated names of the
% 7 genes.  This is useful for use with the "legend" command, among other
% things.
%
% TimePoints is a cell array of strings naming the time points, again,
% mainly meant for use with the "legend" command.

global GeneExpr GeneNames TimePoints

GeneExpr = zeros(7,58,7);

GeneExpr(:,:,1) = [ 0.2115     0.20774     0.20286     0.19599     0.18911     0.18224     0.17437     0.16502     0.15566      0.1463     0.13677     0.12701     0.11726      0.1075    0.098543    0.090442    0.082342    0.074241    0.066692    0.059651     0.05261    0.045569    0.040169    0.036059     0.03195     0.02784    0.024491    0.021649    0.018806    0.015964    0.014124    0.012846    0.011569    0.010291     0.00931   0.0084682   0.0076265   0.0067847   0.0058217   0.0048116   0.0038015   0.0027914   0.0023308   0.0020436   0.0017564   0.0014692   0.0011266  0.00077005  0.00041355  5.7041e-05           0           0           0           0  0.00015686  0.00033512  0.00051337  0.00069162
      0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0
      0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0
      0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0
0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
0.16816      0.1805     0.19538     0.21479      0.2342     0.25361     0.27073     0.28442      0.2981     0.31179     0.32662     0.34291      0.3592     0.37549     0.39096     0.40553     0.42011     0.43469     0.44847     0.46151     0.47455     0.48759     0.49537       0.499     0.50264     0.50627     0.51111     0.51674     0.52238     0.52801     0.52715     0.52263     0.51812      0.5136     0.51024     0.50743     0.50462     0.50181      0.4999     0.49835     0.49679     0.49524     0.49191     0.48803     0.48415     0.48027     0.47669      0.4732      0.4697     0.46621     0.46734     0.46935     0.47136     0.47337     0.46542     0.45613     0.44683     0.43753
      0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0];
  
GeneExpr(:,:,2) = [  0.3565     0.35072     0.34351     0.33377     0.32352     0.31245     0.30065     0.28775     0.27222     0.25305     0.23423     0.21585     0.19575     0.17365     0.15269     0.13296      0.1151    0.099093    0.085541    0.074252    0.064364    0.055669    0.047396    0.039453    0.032936    0.027449     0.02314    0.019614     0.01631    0.013141    0.010606    0.008427   0.0064706   0.0046286   0.0034466   0.0025752   0.0017453  0.00093325   0.0004492  9.2692e-05           0           0           0           0           0           0           0           0           0           0           0           0    0.001039   0.0022472   0.0034031   0.0045518   0.0056827   0.0068116
0.037564     0.04101    0.044563    0.048307    0.052253    0.056531     0.06119    0.066418    0.071872    0.077635    0.083503    0.089504    0.095515     0.10154     0.10608     0.10903     0.11162     0.11384      0.1156     0.11695     0.11577     0.11244     0.10815      0.1031    0.097994    0.092844    0.086055    0.078172    0.071763    0.066257    0.061752    0.057811    0.052876    0.047429     0.04244    0.037667     0.03417    0.031218    0.027897    0.024431    0.022093    0.020152    0.018618    0.017211    0.014955     0.01246    0.011359    0.010606   0.0099671   0.0093531   0.0092216    0.009182   0.0097897    0.010503    0.011425    0.012376    0.012239    0.011981
0.084392    0.079777    0.074556    0.068258    0.061561    0.054213    0.046793    0.039267    0.033089    0.028771    0.025664      0.0241    0.022745    0.021636    0.020745    0.020092    0.020171    0.020983    0.021671    0.022246     0.02297    0.023821    0.023963    0.023547    0.023981    0.025031     0.02658    0.028461    0.029643    0.030396    0.032872    0.036318    0.040549    0.045184    0.046896    0.047232    0.050799    0.055751    0.057351    0.057648    0.062298    0.068478    0.074296        0.08    0.084221    0.088024    0.089133    0.089569     0.09106    0.092783    0.096053    0.099619    0.099607    0.099012    0.099429    0.099984    0.099932    0.099813
0.011907    0.013016    0.012528   0.0092002   0.0082812    0.011292    0.013851    0.015732    0.016616     0.01612    0.016035    0.016471    0.017535    0.019337    0.022813    0.028101    0.032617     0.03636    0.042349     0.05041    0.060021    0.070954    0.084538      0.1002       0.117      0.1346     0.15237     0.17023      0.1857      0.1997     0.20276     0.19965     0.19614     0.19244     0.19155     0.19199     0.17909     0.16047     0.14621     0.13366     0.11812     0.10155    0.091638    0.083834    0.072524    0.060225    0.056339    0.054556    0.053488    0.052577    0.055143    0.058372    0.061174    0.063907    0.062057    0.059581    0.057123    0.054667
 0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
  0.1341     0.14133     0.15092     0.16468     0.17942     0.19576     0.21138     0.22592     0.23977     0.25266     0.26561     0.27862     0.29186     0.30537     0.31986     0.33543     0.34844      0.3589     0.36866     0.37777       0.388     0.39919     0.40823     0.41558     0.42346     0.43172     0.43813     0.44332     0.44446     0.44311     0.44704     0.45393     0.45085     0.44263     0.44291     0.44719     0.44796     0.44722     0.44264     0.43656     0.43379     0.43219     0.43341     0.43553     0.42956      0.4213     0.42551     0.43284     0.43026     0.42551     0.42436     0.42391     0.42469     0.42569     0.42631     0.42688     0.41662     0.40515
       0           0           0           0           0           0           0           0           0           0           0  0.00045469   0.0010143   0.0020636   0.0023317   0.0025882   0.0013058   0.0022734   0.0021102   0.0032528   0.0048384   0.0016089   0.0019237  0.00083943           0           0  0.00030313           0   0.0013408           0  0.00074616  1.1659e-05  0.00072284           0  0.00080445  0.00061791   0.0035792    0.001539   0.0016672   0.0062141   0.0076714   0.0091404    0.013233    0.018188    0.017173    0.023131    0.032038    0.043918    0.057093    0.079722     0.10305     0.12734     0.14649     0.16874     0.18228     0.19338      0.1979     0.20086];

GeneExpr(:,:,3) = [  0.60741     0.59961     0.58812     0.57494     0.55761     0.54439     0.52608     0.49937     0.46439     0.42933     0.38745     0.34361     0.29635     0.25216      0.2049     0.16239     0.12875     0.10106    0.075608    0.059412    0.044863    0.034392    0.026078    0.019647    0.015294    0.011569   0.0081961   0.0057255   0.0037647   0.0026275   0.0011765  0.00082353           0           0           0           0           0           0           0           0           0           0           0           0           0  0.00078431   0.0014118   0.0043529   0.0056078   0.0074902    0.010039    0.012196    0.014588    0.016078    0.017608    0.021608    0.022588     0.02298
 0.037255    0.044392    0.053333    0.061961       0.074    0.089412     0.10965     0.12824     0.15349     0.17992     0.21051     0.24004     0.26204      0.2878     0.31173     0.33388     0.35651     0.36671     0.36855     0.37957       0.378     0.37427     0.36498     0.34592     0.32314       0.294     0.26373     0.22937     0.19584     0.16231     0.13173     0.10471    0.082745    0.066706    0.052902    0.043961    0.035725    0.030275    0.023804    0.019216    0.014863    0.013255     0.01102   0.0093725   0.0078039   0.0072941   0.0062745   0.0032549   0.0030588    0.003098   0.0024314   0.0021569   0.0021961   0.0031373   0.0023922   0.0032157   0.0036471   0.0031765
  0.24071      0.2209     0.19365     0.16545     0.14051     0.10851    0.085529    0.060471    0.043804    0.031373    0.020784    0.013882    0.010392   0.0071373   0.0042353   0.0021569   0.0021961  0.00058824  0.00039216  0.00039216           0           0           0           0           0           0           0   0.0005098   0.0026667   0.0053725    0.011529     0.01651    0.024235    0.033882        0.05    0.064824    0.082039     0.10529     0.12549     0.14894     0.16914     0.18961     0.20718     0.22337     0.23098     0.23851      0.2471     0.24333     0.23898     0.23047     0.21498     0.19988     0.18573     0.17302     0.16063     0.14765     0.13318     0.12165
0.0053333   0.0060784   0.0041961   0.0031373   0.0053333   0.0049412   0.0052941   0.0056078       0.006   0.0051765   0.0058039   0.0059216   0.0064706   0.0075294    0.010196    0.013608    0.013176       0.018    0.024706    0.036353    0.043647     0.06251    0.083882     0.11314     0.14314     0.17537     0.20741     0.23965     0.26035     0.28635     0.29349     0.29678     0.29725     0.29561     0.28686     0.26604     0.24518     0.22224     0.19663     0.16675       0.138     0.11396    0.091451    0.072118    0.056118    0.047725    0.034471    0.029765    0.023765    0.019608    0.018314    0.014941    0.013373    0.013216    0.013059    0.011686    0.012941    0.011255
  0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
  0.12467     0.14059     0.15373     0.16412     0.18082     0.20243     0.21776     0.22882     0.25078     0.26769     0.28808     0.29875     0.31455      0.3329     0.34898     0.37753     0.38753     0.39776     0.41718     0.43941     0.44584     0.45871      0.4829     0.48647     0.49592     0.50537     0.51894     0.52478     0.53584     0.53867     0.54463     0.54882     0.55008     0.55396     0.56337     0.56314     0.56553     0.56471     0.55855      0.5658     0.56533     0.56639     0.55584     0.56039     0.54282     0.54141     0.55047     0.53835     0.52933     0.53078     0.52541     0.53043     0.52161     0.52424     0.51871     0.51129     0.52129     0.50169
        0           0           0           0           0           0           0           0           0           0           0  0.00099205    0.002213   0.0045024   0.0050874   0.0056471    0.002849   0.0049603   0.0046041    0.007097    0.010556   0.0035103   0.0041971   0.0018315           0           0  0.00066137           0   0.0029253           0    0.001628  2.5437e-05   0.0015771           0   0.0017552   0.0013482   0.0078092   0.0033577   0.0036375    0.013558    0.016738    0.019943    0.028871    0.039682    0.037469    0.050467    0.069901    0.095822     0.12457     0.17394     0.22484     0.27783     0.31962     0.36815     0.39771     0.42193     0.43177     0.43823];

GeneExpr(:,:,4) = [ 0.72239       0.716     0.70416     0.68424     0.66533     0.64412     0.60012     0.55851     0.50965     0.44741     0.38878     0.32184     0.25906     0.20212     0.15729      0.1182    0.092118    0.070824    0.053294    0.041059    0.031529    0.022549    0.017216    0.013412   0.0086667   0.0062745   0.0036863   0.0030196   0.0014118           0           0           0           0           0           0           0           0           0           0           0  0.00039216   0.0017255   0.0038431   0.0064314    0.011137    0.016745    0.021647    0.028392    0.035373    0.039765    0.046078    0.051647    0.057373    0.061412    0.065373    0.068902    0.070353    0.064941
0.034196    0.042392    0.055255    0.070314    0.092745     0.12055     0.15443     0.18569     0.23827     0.28341     0.33733      0.3829     0.42553     0.46063     0.49133     0.51247     0.51957     0.52533     0.53075     0.52302     0.50651     0.49094     0.46302     0.42612     0.38502     0.33722     0.28137     0.23341     0.18506     0.14537     0.11024    0.081098    0.059882    0.044275    0.030941     0.02298    0.015647    0.010627   0.0068235   0.0040392   0.0017647  0.00015686           0           0           0           0           0           0           0           0           0  0.00015686   0.0010196   0.0013725   0.0024314   0.0023529   0.0041569   0.0041176
 0.36169     0.31353     0.26114     0.20745     0.15361     0.11118    0.078235    0.055373     0.03651    0.025176       0.018    0.012941   0.0091373   0.0058039   0.0045882    0.003098   0.0026667  0.00082353  0.00062745   0.0011765   0.0012941  0.00090196  0.00039216  0.00039216  0.00031373  0.00090196   0.0017647   0.0023922   0.0052157   0.0081569    0.014627    0.024588    0.040196    0.058902    0.084902     0.12161     0.15655     0.19459     0.22882     0.26094       0.284     0.30608     0.32227     0.32663       0.328     0.31349     0.29294     0.27443     0.24929     0.22102     0.20573       0.184     0.16302     0.14682     0.13098     0.11729     0.10216    0.092353
       0           0           0           0           0           0           0           0           0           0           0           0           0  0.00035294   0.0005098   0.0035686   0.0088627    0.013059    0.020549    0.039922    0.065294    0.099412     0.14404     0.19827     0.25133     0.31098     0.36914     0.40918     0.45118     0.45875     0.46686     0.46373     0.44659     0.40537     0.36596     0.32051     0.27031     0.21694     0.17671     0.14455        0.11    0.082706    0.064745    0.048706    0.036902     0.02902     0.02302    0.019137    0.013333    0.012078    0.010784    0.009451   0.0063922   0.0061961   0.0072157   0.0050588   0.0043529   0.0033333
 0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
  0.1191     0.13561     0.13612     0.15494     0.16357     0.18412     0.19416     0.21302     0.22584     0.23271     0.24776     0.26859     0.27965     0.30345     0.31243     0.33404     0.34173     0.34557     0.36851      0.3831     0.39376     0.40239     0.41863     0.42447     0.44447     0.45667     0.45918     0.47325     0.46839     0.47416     0.47251     0.49129     0.49576     0.50314     0.50059     0.50302     0.49965     0.50573      0.5109     0.50651     0.51639     0.50498     0.51345     0.51404     0.51196      0.5091     0.52318      0.5171     0.52067     0.52204      0.5178     0.52247     0.51706     0.51722     0.51298     0.50216     0.50788     0.49192
       0           0           0           0           0           0           0           0           0           0           0   0.0012401   0.0027663    0.005628   0.0063593   0.0070588   0.0035612   0.0062003   0.0057552   0.0088712    0.013196   0.0043879   0.0052464   0.0022893           0           0  0.00082671           0   0.0036566           0    0.002035  3.1797e-05   0.0019714           0    0.002194   0.0016852   0.0097615   0.0041971   0.0045469    0.016948    0.020922    0.024928    0.036089    0.049603    0.046836    0.063084    0.087377     0.11978     0.15571     0.21742     0.28105     0.34728     0.39952     0.46019     0.49714     0.52741     0.53971     0.54779];

GeneExpr(:,:,5) = [  0.80533      0.7989      0.7878     0.77976     0.76243     0.74714     0.72796     0.69169     0.64169     0.57063     0.49596     0.40114     0.31565     0.24086     0.17886     0.13298    0.097961    0.070824    0.052471    0.037922    0.027098    0.019451       0.014    0.010235   0.0072941   0.0050588   0.0028235   0.0013725  0.00047059           0           0           0           0           0           0           0           0           0           0  0.00070588   0.0034902   0.0076863    0.014157    0.024824    0.039373    0.056784    0.076353    0.091922     0.10733     0.11902     0.12918     0.13541     0.13675     0.13604     0.13161     0.12447     0.11322     0.10482
 0.020863       0.026    0.035137    0.046039    0.064275     0.08898     0.12396     0.17133     0.22318     0.29161     0.36392     0.43396     0.49961     0.55239      0.5858     0.60631     0.61812     0.62282     0.61969     0.60761     0.59718     0.56976     0.53553      0.4822     0.41867     0.34533     0.27541     0.20455     0.14765      0.1069    0.074588    0.053412    0.038314    0.026824    0.018314    0.013569   0.0090196   0.0063922   0.0032157       0.002  0.00086275  0.00082353           0           0           0           0  0.00031373   0.0007451    0.001098  0.00098039   0.0016471   0.0018824   0.0016471   0.0026275   0.0022353   0.0021961   0.0023922   0.0023922
     0.54     0.48486     0.41973     0.33694     0.25514     0.18035     0.12114    0.079686    0.055059    0.036353    0.021843    0.015451    0.010667   0.0077647   0.0065882   0.0045882   0.0030588   0.0020784   0.0019216       0.002   0.0017647  0.00070588  0.00035294           0  0.00039216           0  0.00039216   0.0020784    0.004902    0.012078    0.024392    0.044235    0.078588     0.12208     0.16882      0.2291     0.28345     0.33275     0.37388     0.39373     0.40396     0.40502     0.38953     0.36427     0.33075     0.29408     0.25957     0.22741     0.19569     0.16643     0.14812     0.12675     0.11129    0.092275    0.081333    0.069176    0.058392     0.05102
0.0023529   0.0016471           0           0           0  0.00019608           0    0.001451   0.0018431   0.0012157           0  0.00094118           0           0   0.0011765   0.0040392   0.0075294       0.014    0.024039    0.038706    0.062549     0.10604     0.15808     0.22267     0.29827     0.38573       0.468     0.53624     0.59031     0.61396     0.62031     0.60812     0.57776      0.5282     0.45965     0.38451     0.30745     0.23298     0.17376     0.13114     0.10161    0.076706    0.057961    0.040118    0.032824    0.024314    0.020275     0.01702     0.01451    0.011451     0.01098   0.0089804   0.0088235    0.007098   0.0066667   0.0061176   0.0056471   0.0045882
  0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
 0.051137    0.066745    0.076863    0.086275    0.099608     0.11024       0.124     0.14031     0.14918     0.16498     0.17839     0.19227      0.2138     0.23031      0.2329     0.24678     0.26173     0.27659     0.29271         0.3     0.31035     0.32525     0.33808     0.33729     0.34271     0.35467     0.36051     0.36875     0.37004     0.37522     0.37235     0.38455     0.39329     0.39663     0.39718     0.40604     0.40455     0.41012     0.41847     0.42004     0.42651     0.44192     0.44427     0.44502     0.43851     0.44937     0.45004      0.4471     0.45408     0.44443     0.45098     0.45851       0.458     0.44965     0.42941     0.42365       0.406     0.39282
        0           0           0           0           0           0           0           0           0           0           0   0.0015294   0.0034118   0.0069412   0.0078431   0.0087059   0.0043922   0.0076471    0.007098    0.010941    0.016275   0.0054118   0.0064706   0.0028235           0           0   0.0010196           0   0.0045098           0   0.0025098  3.9216e-05   0.0024314           0   0.0027059   0.0020784    0.012039   0.0051765   0.0056078    0.020902    0.025804    0.030745     0.04451    0.061176    0.057765    0.077804     0.10776     0.14773     0.19204     0.26816     0.34663     0.42831     0.49275     0.56757     0.61314     0.65047     0.66565     0.67561];

GeneExpr(:,:,6) = [  0.81525      0.8109     0.81255     0.80961     0.80635     0.80796       0.794     0.76839     0.72973     0.65424     0.57286     0.47467     0.37447     0.27937     0.20976     0.15165     0.10922    0.075882    0.052275    0.037059    0.026706    0.019529     0.01451    0.010902   0.0070196    0.005451   0.0034118   0.0018431  0.00031373           0           0           0           0           0           0           0           0           0  0.00070588   0.0054902    0.013922    0.028196    0.051608    0.081255      0.1162     0.14984     0.18435     0.21902      0.2442     0.25541     0.25882     0.25941     0.24816     0.23518     0.20965     0.18702     0.15941     0.14008
 0.026392    0.033098    0.045059    0.060824    0.086314     0.12369     0.17529     0.25008     0.33176     0.44043     0.52776     0.60494      0.6638     0.70027     0.72439      0.7331     0.73596     0.73827     0.73243      0.7138     0.68671     0.64878      0.5869     0.51792     0.42894     0.32973     0.23937     0.16235     0.11388    0.078118    0.052196    0.035961    0.025725     0.01898    0.012745    0.008549   0.0064706   0.0048627   0.0039216   0.0021961   0.0011373  0.00062745  0.00019608           0  0.00062745  0.00082353  0.00070588   0.0014118   0.0012941   0.0017255    0.001451   0.0013333  0.00098039  0.00054902  0.00043137  0.00019608  0.00035294           0
  0.60616     0.52102     0.40902     0.30529      0.2022     0.13416    0.080275    0.049137    0.030549    0.021412    0.014118   0.0094118   0.0058039   0.0053333   0.0031373    0.003098   0.0016471  0.00058824           0           0           0           0           0           0           0           0           0   0.0018431   0.0093333     0.02451    0.056314     0.10196     0.17129     0.25808     0.34439     0.42984     0.49475     0.54318     0.56796     0.58231     0.57561     0.54549     0.49576     0.42878     0.35961     0.29325     0.23796     0.18337     0.15122     0.12827     0.10424    0.088588    0.073882    0.061255    0.049961    0.040275    0.034667    0.026902
        0   0.0012157           0           0           0           0   0.0019216   0.0031765   0.0036471    0.001451   0.0023922           0  0.00078431   0.0039608       0.008    0.012941    0.019098    0.025882    0.043843    0.073451     0.11929     0.18373     0.26886     0.37169     0.49227     0.60671     0.68918     0.74078     0.76424     0.76016     0.73412     0.68537     0.61627     0.52545     0.42949     0.33345     0.25216     0.19122     0.14322     0.10576    0.076902    0.059922    0.047451    0.036431    0.030078    0.027059    0.024353    0.023098    0.021373    0.018549    0.018118    0.014078    0.014627    0.011725    0.010471   0.0090588    0.007098   0.0063529
  0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
 0.060118    0.071373    0.083725    0.093843     0.10502     0.11486     0.12776     0.14341     0.15118     0.17196     0.18075     0.19898     0.21718     0.22541      0.2331     0.23086     0.24259     0.25306     0.26839      0.2718     0.27675     0.28349     0.29165     0.29635     0.30329     0.30525     0.30902     0.31016     0.31647      0.3111     0.31224     0.32161     0.32643     0.33608     0.33325     0.32961        0.33     0.33796     0.34008     0.34514       0.348     0.35231       0.358     0.36988     0.37831     0.37263     0.37812     0.37671      0.3871     0.37941     0.38655     0.38875     0.38886     0.37133      0.3669     0.34271     0.33184     0.31455
0.0066275   0.0039216   0.0088627    0.006902       0.008  0.00035294  0.00058824   0.0027843           0           0   0.0011765           0           0   0.0016863           0    0.001451           0   0.0039608   0.0088627    0.015176   0.0018431           0   0.0020392           0   0.0010588   0.0035686   0.0086667   0.0043529   0.0026667           0           0           0           0           0           0           0           0           0           0   0.0021961   0.0095686    0.012627    0.023373    0.025176     0.03451    0.043843     0.05949    0.076863     0.10882     0.13992     0.17929     0.24435     0.30643     0.38957     0.43569     0.50145     0.50376     0.55843];

GeneExpr(:,:,7) = [ 0.80996     0.80506     0.79922     0.80267     0.80573     0.80369     0.79533     0.78302       0.758      0.7042     0.61706     0.50588     0.39565     0.27871     0.18788     0.12867       0.088    0.060314    0.042157    0.030667    0.021098    0.013647    0.010431   0.0072549   0.0056471       0.004   0.0031373   0.0020392   0.0019608   0.0011373  0.00058824           0           0           0           0           0           0           0           0   0.0036078     0.01098    0.026157    0.047882    0.082745     0.12071     0.16196     0.19918     0.23361      0.2602     0.28549         0.3     0.30788     0.30373     0.28851     0.26071     0.23125     0.20004     0.16573
0.015804    0.018784    0.024549    0.032314    0.045294    0.067333    0.097961     0.14494     0.21671     0.31863     0.45776     0.57275     0.65882     0.72839     0.76204     0.77416     0.79173     0.80106     0.80031      0.7902     0.77063     0.73592     0.68282     0.60949     0.52882     0.42145     0.30729     0.21455     0.14502    0.092235    0.063647    0.042039    0.028471    0.019804    0.014235   0.0098824   0.0072157   0.0064706   0.0047059   0.0036471   0.0020784   0.0015294   0.0012941   0.0012549   0.0017647   0.0018824   0.0023137   0.0024314   0.0028235   0.0029412   0.0029804    0.002549       0.002   0.0014118   0.0011373   0.0005098           0           0
 0.80035     0.73475     0.63675     0.50769     0.35525     0.20596     0.11616    0.062784    0.031451    0.018863    0.011059   0.0080392   0.0064706   0.0048235   0.0053725   0.0048627   0.0041176    0.003098   0.0031765   0.0027059   0.0017647   0.0015686  0.00027451           0           0           0           0           0  0.00086275   0.0098824    0.029882    0.072039     0.14843     0.26902     0.38243     0.49863     0.57729     0.63286     0.65847     0.66031     0.64118     0.59651     0.52867     0.44741     0.35078     0.26075      0.1898     0.13957     0.11137    0.087137    0.070118    0.056824     0.04902    0.038941    0.031216    0.026157     0.02102    0.015333
       0           0           0           0           0  0.00082353    0.003098   0.0032941   0.0033725   0.0014118           0           0           0   0.0028627   0.0077255    0.011686    0.013765    0.019529    0.033922    0.062549     0.10902       0.166     0.24569     0.34788     0.48624     0.61875     0.71584     0.77839     0.79957     0.81337     0.79761     0.76871     0.70447     0.60824     0.48314     0.36396     0.26302       0.192     0.13898      0.1069    0.077137    0.059412    0.047059    0.037294    0.030706     0.02698    0.024549    0.022196    0.019765    0.017255    0.013804    0.012353   0.0099216   0.0081961   0.0058039   0.0050196   0.0047843   0.0021961
 0.17518     0.16725     0.15792     0.15027     0.14267     0.13459     0.12745     0.12251     0.11478     0.10859     0.10286    0.097843    0.092196    0.087176    0.081608    0.076706     0.07298    0.068745    0.065451    0.061255    0.057294    0.054157    0.051255    0.048235    0.045451    0.042549    0.040353    0.037843    0.035137    0.033255    0.031647    0.029961    0.028118    0.026706    0.024941    0.023333    0.022431    0.021686    0.020314    0.019294    0.017882     0.01698    0.016039    0.014941     0.01451    0.013608    0.012706    0.012314    0.011765    0.011333     0.01102    0.010784    0.010588    0.010078   0.0098824   0.0096078   0.0090588   0.0084314
   0.018    0.029647    0.035686    0.052902    0.051176    0.060314    0.074118    0.078314    0.095137     0.11306     0.12737     0.13922     0.14455     0.15988     0.15847     0.16122     0.17227     0.18098     0.19839     0.20741     0.22353     0.22576     0.20886     0.22447     0.23694     0.24165     0.25094     0.25537     0.25553     0.25576     0.25643      0.2609      0.2669     0.27125      0.2798     0.28639     0.28325     0.27961     0.28553     0.30722     0.32408     0.32918     0.34455     0.34545     0.34388      0.3451      0.3482     0.35047     0.35682     0.36667     0.37157     0.36459     0.37541     0.36318     0.34192     0.32267     0.31247     0.30063
       0           0           0           0   0.0021176           0   0.0022353           0           0           0           0           0           0           0           0           0           0           0           0    0.012824   0.0013333           0           0   0.0047843           0   0.0056863    0.010588   0.0090196   0.0047059   0.0034902   0.0027843   0.0051373   0.0043529   0.0076471    0.012745   0.0060784        0.01   0.0062353   0.0096471   0.0055294    0.008902    0.017608    0.026667    0.033529    0.033137    0.045176    0.059765    0.085922     0.12675     0.17098     0.21627     0.29647     0.36894     0.45514     0.52604     0.58243     0.63945     0.66122];

% Normalizing GeneExpr:
for g=1:7
    GeneExpr(g,:,:) = GeneExpr(g,:,:)/max(max(max(GeneExpr(g,:,:))));
end
   
GeneNames = {'Hb','Kr','Gt','Kni','Bcd','Cad','Tll'};

TimePoints = {'t=0','t=1','t=2','t=3','t=4','t=5','t=6','t=7'};

