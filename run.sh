#rootcint -f dict.C RooDoubleCBFast.h
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  makePlots.C -o makePlots
g++ -g -O3 `root-config --libs --cflags` -lGenVector -lRooFit -lRooFitCore AtlasLabels.C  AtlasUtils.C  RooDoubleCBFast.C  dict.C  treatPlots.C -o treatPlots
./makePlots
./treatPlots
mv *.png *.pdf ~/eos/www/low_pt_jets/
