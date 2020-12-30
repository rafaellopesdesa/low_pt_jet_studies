#rootcint -f dict.C RooDoubleCBFast.h
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makePlots.C -o makePlots
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  treatPlots.C -o treatPlots
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makeCutflow.C -o makeCutflow
./makePlots
./treatPlots
./makeCutflow
rm -rf ~/eos/www/low_pt_jets/*.png 
rm -rf ~/eos/www/low_pt_jets/*.pdf 
mv *.png *.pdf ~/eos/www/low_pt_jets/
