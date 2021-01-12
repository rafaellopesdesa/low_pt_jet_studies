#rootcint -f dict.C RooDoubleCBFast.h
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makePlots.C -o makePlots
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  treatPlots.C -o treatPlots
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makeCutflow.C -o makeCutflow
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makePlotsZee.C -o makePlotsZee
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makePlotsData18.C -o makePlotsData18
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  treatPlotsZee.C -o treatPlotsZee
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  treatPlotsData18.C -o treatPlotsData18
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  compareDataMC.C -o compareDataMC
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore -lMinuit AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  fitScale.C -o fitScale
./makePlots
./treatPlots
./makeCutflow
rm -rf ~/eos/www/low_pt_jets/*.png 
rm -rf ~/eos/www/low_pt_jets/*.pdf 
mv *.png *.pdf ~/eos/www/low_pt_jets/
