void generate_lymap(){
double x_start=-6;
double x_end=6;
const int nbins = 50;
double step;

step = (x_end-x_start)/nbins;

double xx[nbins];
double ly[nbins];


gRandom->SetSeed(0);

for(int i=0;i<nbins;i++){
	xx[i] = x_start+i*step;
	ly[i] = gRandom->Gaus(180,0.1*180); 
}

for(int i=0;i<nbins;i++){
	cout<<xx[i]<<", ";
}
cout<<endl;
for(int i=0;i<nbins;i++){
	cout<<ly[i]<<", ";
}
cout<<endl;

}
