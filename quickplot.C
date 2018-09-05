#include <float.h>

TTree *text2tree(const char *filename){
    /*util function to contruct and fill a TTree from a text file.*/



}

TH1D *quickplot_file(const char *filename, int xcol=0, int ycol=1, int skiplines=0, int bins=0, string xfunc="", string yfunc=""){
    FILE *fp;
    fp=fopen(filename,"r");
    char line[256];
    std::vector<double> X,Y;
    for (int ii=0;ii<skiplines;ii++){
        fgets(line,255,fp);
    }
    fgets(line,255,fp);
    double xmax=-DBL_MAX,xmin=DBL_MAX;
    double ymax=-DBL_MAX,ymin=DBL_MAX;
    int linec=0;
    int I=0;
    while(!feof(fp)){
        std::vector<double> buffer;
        vector<double> tmp;
        if(line[0]=='#' || strlen(line)<3) {
            fgets(line,255,fp);       
            continue;
        }
        string tmpstr;
        stringstream iss(line);
        while (std::getline(iss,tmpstr,' ')){
            if (!tmpstr.empty()){
                tmp.push_back(atof(tmpstr.c_str()));
            }
        }
        if(xcol>=0){
            X.push_back(tmp[xcol]);
            if (xmax<tmp[xcol])xmax=tmp[xcol];
            if (xmin>tmp[xcol])xmin=tmp[xcol];
        }else{
            X.push_back(I++);
            xmin=0;
            xmax=I;
        }
        Y.push_back(tmp[ycol]);
       
        if (ymax<tmp[ycol])ymax=tmp[ycol];
        if (ymin>tmp[ycol])ymin=tmp[ycol];

        linec++;
        fgets(line,255,fp);
    }
    
    fclose(fp);
    if (bins) linec=bins;
    TH1D *h1 = new TH1D("yxhist", "Y=f(X)", linec, xmin,xmax);
    cout << "lines " << linec << " xmx " << xmax << " xmn " << xmin <<endl;
    for (int ii=0;ii<X.size();ii++){
        h1->Fill(X[ii],Y[ii]);
    }
    h1->Draw();
    return h1;
}

TH2D *quicksplot_file(const char *filename, int xcol=0, int ycol=1, int zcol=2, int skiplines=0, int bins1=0, int bins2=0){
    FILE *fp;
    fp=fopen(filename,"r");
    char line[256];
    std::vector<double> X,Y,Z;
    for (int ii=0;ii<skiplines;ii++){
        fgets(line,255,fp);
    }
    fgets(line,255,fp);
    double xmax=-DBL_MAX,xmin=DBL_MAX;
    double ymax=-DBL_MAX,ymin=DBL_MAX;
    double zmax=-DBL_MAX,zmin=DBL_MAX;
    int linec=0;
    while(!feof(fp)){
        std::vector<double> buffer;
        vector<double> tmp;
        if(line[0]=='#' || strlen(line)<3) {
            fgets(line,255,fp);       
            continue;
        }
        string tmpstr;
        stringstream iss(line);
        while (std::getline(iss,tmpstr,' ')){
            if (!tmpstr.empty()){
                tmp.push_back(atof(tmpstr.c_str()));
            }
        }
        X.push_back(tmp[xcol]);
        Y.push_back(tmp[ycol]);
        Z.push_back(tmp[zcol]);

        if (xmax<tmp[xcol])xmax=tmp[xcol];
        if (ymax<tmp[ycol])ymax=tmp[ycol];
        if (zmax<tmp[zcol])zmax=tmp[zcol];

        if (xmin>tmp[xcol])xmin=tmp[xcol];
        if (ymin>tmp[ycol])ymin=tmp[ycol];
        if (zmin>tmp[zcol])zmin=tmp[zcol];

        linec++;
        fgets(line,255,fp);
    }
    
    fclose(fp);
    TH2D *h2 = new TH2D("yxzhist", "Z=f(X,Y)", bins1, xmin,xmax, bins2, ymin,ymax);
    for (int ii=0;ii<X.size();ii++){
        h2->Fill(X[ii],Y[ii],Z[ii]);
    }
    h2->SetStats(kFALSE);
    h2->Draw("colz");
    return h2;
}


TH2D *readstuff(const char *filename, double emin, double emax){
    
    TH2D *h2 = new TH2D("spatres_40m", "Spatially resolved intensity", 41, -8e-3,8e-3, 41, -2.1e-3,2.1e-3);
    h2->SetXTitle("x / m");
    h2->SetYTitle("y / m");
    h2->SetZTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 

    FILE *fp;
    fp=fopen(filename,"r");
    char line[256];
    fgets(line,255,fp);
    while(!feof(fp)){
        if(line[0]=='#' || strlen(line)<3) {
            fgets(line,255,fp);       
            continue;
        }
        double energy,xi,yi,II;
        int li;
        sscanf(line,"%lf %lf %lf %lf",&energy,&yi,&xi,&II);
        if(energy<emax && energy>=emin){
            //cout << II <<endl;
            cout << energy << endl;
            //li=h2->GetBin(xi,yi);
            h2->Fill(xi/40*0.016-0.008,yi/40*4.2e-3-2.1e-3,II);
        }        
        fgets(line,255,fp);       
    }
    fclose(fp);
    return h2;
}

TH1F *readstuff1(const char *filename, double emin, double emax, string axis="x", int xii=0, int yii=0){
    TH1F *h1; 
    Double_t bwx=16.0e-3/40.0;
    Double_t bwy=4.2e-3/40.0;
    //cout <<bwx << " " << bwy <<endl;
    if(axis=="x"){ 
        h1 = new TH1F("@40m_x", "X-ray intensity, cut along x", 41, -41.0*bwx*0.5,41*bwx*0.5);
        h1->SetXTitle("x / m");
        h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
    }else{
        h1 = new TH1F("@40m_y", "X-ray intensity, cut along y", 41, -41*bwy*0.5,41*bwy*0.5);
        h1->SetXTitle("y / m");
        h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
    }
    FILE *fp;
    fp=fopen(filename,"r");
    char line[256];
    fgets(line,255,fp);
    while(!feof(fp)){
        if(line[0]=='#' || strlen(line)<3) {
            fgets(line,255,fp);       
            continue;
        }
        double energy,II;
        int li,;
        int xi,yi;
        li=sscanf(line,"%lf %ld %ld %lf",&energy,&yi,&xi,&II);
        if(energy<=emax && energy>=emin){
            if ( ((xii && xi==xii)|| (yii && yi==yii) ) && ((xi==xii || yi==yii)) ) cout << "choosing line number xi "<<xi <<" and yi "<<yi<<endl;
            if (axis=="x" && (yii!=0 && yi==yii)){
                h1->Fill(xi*0.016/40-0.008,II);
                cout << xi*0.0016/40-0.008 << " " << II<<endl;
                //if ((xii || yii) && (xi==xii || yi==yii)) cout << "choosing line number xi "<<xi <<" and yi "<<yi<<endl;
                //h1->Fill(xi,II);
            }else if (axis=="y" && (xii!=0 && xi==xii)){
                h1->Fill(yi*4.2e-3/40-2.1e-3,II);
                cout << yi*4.2e-3/40-2.1e-3 << " " << II<<endl;
            }

        }        
        fgets(line,255,fp);       
    }
    fclose(fp);
    return h1;
}

TH1F *readstuff1(const char *filestem, string axis="x", int skiplines=10,int istart=0, int iend=10){
    TH1F *h1;
    double xmax,xmin,ymax,ymin;
    for (int ii=istart;ii<iend; ii++){ 
        std::vector<double> X,Y;
        int linec=0;
        const char filename[1024];
        snprintf(filename,1024,"%s%d.dt%s",filestem,ii,axis.c_str());
        cout << filename <<endl;
        FILE *fp;
        fp=fopen(filename,"r");
        char line[256];
        fgets(line,255,fp);
        for (int jj=0;jj<skiplines;jj++){
            fgets(line,255,fp);
        }
        while(!feof(fp)){
            if(line[0]=='#' || strlen(line)<3) {
                fgets(line,255,fp);       
                continue;
            }
            vector<double> tmp;
            string tmpstr;
            stringstream iss(line);
            while (std::getline(iss,tmpstr,' ')){
                if (!tmpstr.empty()){
                    tmp.push_back(atof(tmpstr.c_str()));
                }
            }
            X.push_back(tmp[0]);
            Y.push_back(tmp[1]);

            if (xmax<tmp[0])xmax=tmp[0];
            if (ymax<tmp[1])ymax=tmp[1];

            if (xmin>tmp[0])xmin=tmp[0];
            if (ymin>tmp[1])ymin=tmp[1];

            linec++;
            fgets(line,255,fp);
        }
        //cout << "lines " << linec << "xmx " << xmax << " xmn " << xmin <<endl;
        //cout << "lines " << linec << "ymx " << ymax << " ymn " << ymin <<endl;
        if(ii==istart){
            if(axis=="x"){ 
                h1 = new TH1F("@40m_x", "X-ray intensity, cut along x", linec,xmin,xmax);
                h1->SetXTitle("x / m");
                h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
                h1->Draw();
            }else{
                h1 = new TH1F("@40m_y", "X-ray intensity, cut along y", linec, xmin,xmax);
                h1->SetXTitle("y / m");
                h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
            }
        }
        for (int jj=0;jj<X.size();jj++){
            h1->Fill(X[jj],Y[jj]);
        }
    }
    h1->SetDrawOption("nostats");
    h1->Draw();
    return h1;
}

TH1F *readstuff1(const char *filestem, string axis="x", int skiplines=10,int istart=0, int iend=10, double limit1, double limit2,int bins=10){
    TH1F *h1;
    double xmax,xmin,ymax,ymin;
    for (int ii=istart;ii<iend; ii++){ 
        std::vector<double> X,Y,Z;
        int linec=0;
        const char filename[1024];
        snprintf(filename,1024,"%s%d.dta",filestem,ii);
        cout << filename <<endl;
        FILE *fp;
        fp=fopen(filename,"r");
        char line[256];
        fgets(line,255,fp);
        for (int jj=0;jj<skiplines;jj++){
            fgets(line,255,fp);
        }
        while(!feof(fp)){
            if(line[0]=='#' || strlen(line)<3) {
                fgets(line,255,fp);       
                continue;
            }
            vector<double> tmp;
            string tmpstr;
            stringstream iss(line);
            while (std::getline(iss,tmpstr,' ')){
                if (!tmpstr.empty()){
                    tmp.push_back(atof(tmpstr.c_str()));
                }
            }
            X.push_back(tmp[0]);
            Y.push_back(tmp[1]);
            Z.push_back(tmp[2]);

            if (xmax<tmp[0])xmax=tmp[0];
            if (ymax<tmp[1])ymax=tmp[1];

            if (xmin>tmp[0])xmin=tmp[0];
            if (ymin>tmp[1])ymin=tmp[1];

            linec++;
            fgets(line,255,fp);
        }
        //cout << "lines " << linec << "xmx " << xmax << " xmn " << xmin <<endl;
        //cout << "lines " << linec << "ymx " << ymax << " ymn " << ymin <<endl;
        if (bins) linec=bins;
        if(ii==istart){
            if(axis=="x"){ 
                h1 = new TH1F("@40m_x", "X-ray intensity, cut along x", linec,xmin,xmax);
                h1->SetXTitle("x / mm");
                h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
            }else{
                h1 = new TH1F("@40m_y", "X-ray intensity, cut along y", linec, ymin,ymax);
                h1->SetXTitle("y / mm");
                h1->SetYTitle("I / pht / s / mm^2 / mrad^2 / 0.1% BW"); 
            }
        }
        for (int jj=0;jj<X.size();jj++){
            if (axis=="x"){
                if (Y[jj]>=limit1 && Y[jj]<limit2) h1->Fill(X[jj],Z[jj]);
            }else{
                if (X[jj]>=limit1 && X[jj]<limit2) h1->Fill(Y[jj],Z[jj]);
            }
        }
    }
    h1->SetStats(kFALSE);//SetDrawOption("nostats");
    h1->Draw();
    return h1;
}

Double_t find_FWHM(TH1F *inh, int doublesided=0){
    double ymax_2=inh->GetMaximum()/2.0;
    int ii=1;
    double ymax,HWHM;
    cout << ymax << " " << inh->GetSize() << endl;
    for(int ii=2;ii<inh->GetSize()-1;ii++){
        double valm1=inh->GetBinContent(ii-1);
        double val=inh->GetBinContent(ii);
        cout <<val <<" " <<valm1<<endl;
        if (val<ymax_2 && valm1>ymax_2){
            double alpha=fabs((ymax_2-valm1)/(val-valm1));
            cout << "alpha= " << alpha << " X[i-1]= " << inh->GetBinCenter(ii-1) << " X[i]= " << inh->GetBinCenter(ii) << endl;
            HWHM=(1-alpha)*inh->GetBinCenter(ii-1) + (alpha)*inh->GetBinCenter(ii);
            cout << "alpha= " << alpha << " HwHM = " << HWHM << " FWHM= " << 2*HWHM << endl; 
            break;   
        }
    }
    return 2*HWHM;
}






Double_t find_FWHM_1dspectra(int skiplines, const char* filename, int ncols=6){
    FILE *fp;
    fp=fopen(filename,"r");
    char line[256];
    std::vector<double> FD,X;
    int li;

    for (int ii=0;ii<skiplines;ii++){
        fgets(line,255,fp);
    }
    while(!feof(fp)){
        if(line[0]=='#' || strlen(line)<3) {
            fgets(line,255,fp);       
            continue;
        }
        double xx,fd,pl[4];
        li=sscanf(line,"%lf %lf %lf %lf %lf %lf",&xx,&fd,&pl,pl+1,pl+2,pl+3);
        X.push_back(xx);
        FD.push_back(fd);
        fgets(line,255,fp);
    }
    double fdmax=-10;
    for(int ii=0;ii<FD.size();ii++){
        if (FD[ii]>fdmax){
            fdmax=FD[ii];
        }
    }
    for(int ii=1;ii<FD.size();ii++){
        if (FD[ii]<fdmax/2.0 && FD[ii-1]>fdmax/2.0){
            double alpha=fabs((fdmax/2.0-FD[ii-1])/(FD[ii]-FD[ii-1]));
            cout << "alpha= " << alpha << " X[i-1]= " << X[ii-1] << " X[i]= " << X[ii] << endl;
            double HWHM=(1-alpha)*X[ii-1] + (alpha)*X[ii];
            cout << "alpha= " << alpha << " HwHM = " << HWHM << " FWHM= " << 2*HWHM << endl; 
        }
    }

}



Double_t find_fwhm(TH1F *h){
    Double_t mm=h->GetMaximum();
    Double_t x1,x2,alpha;

    for (int i=0;i<h->GetSize()-1;i++){
        Double_t val=h->GetBinContent(i);
        Double_t valp1=h->GetBinContent(i+1);
        if (val<mm/2.0 && valp1>mm/2.0){
            alpha=(mm/2.0-val)/(valp1-val);
            x1=(1-alpha)*h->GetBinCenter(i) + (alpha)*h->GetBinCenter(i+1);
            cout << "rising edge in bin " <<i<<", alpha="<<alpha<< ", x1="<<x1<<endl;
        }
       if (val>mm/2.0 && valp1<mm/2.0){
            alpha=(mm/2.0-val)/(valp1-val);
            x2=(1-alpha)*h->GetBinCenter(i) + (alpha)*h->GetBinCenter(i+1);
            cout << "falling edge in bin " <<i<<", alpha="<<alpha<< ", x2="<<x2<<endl;
        }
 
    }
    cout << "FWHM = "<<(x2-x1)<<endl;
    return x2-x1;
}

void runseries (int xi=0, int yi=0){
    TH1F *t1y,*t1x;
    TH2D *t2;

    t1x = readstuff1("rt16_exyI_30kev_40m.dat",30000,30005,"x",0,yi);
    find_fwhm(t1x);
    t1x->Draw();
    t1y = readstuff1("rt16_exyI_30kev_40m.dat",30000,30005,"y",xi,0);
    find_fwhm(t1y);
    //t2 = readstuff1("rt16_exyI_30kev_40m.dat",30000,30005);

    t1x = readstuff1("rt16_exyI_10kev_40m.dat",10000,10005,"x",0,yi);
    find_fwhm(t1x);
    t1y = readstuff1("rt16_exyI_10kev_40m.dat",10000,10005,"y",xi,0);
    find_fwhm(t1y);
    //t2 = readstuff1("rt16_exyI_10kev_40m.dat",10000,10005);
}

TH1F *plotPeak(const char *filestem, double energy=30.0){
    vector<double> xx1,yy1,xx2,yy2;

    for (int ii=0;ii<201;ii++){
        TH1F *hx;
        hx=readstuff1(filestem,"x",10,ii,ii+1,0,8,201);
        xx1=ii*energy*0.1/200.0 + (energy*0.95);
        yy1.push_back(hx->GetSum());
        hy=readstuff1(filestem,"y",10,ii,ii+1,0,8,101);
        xx2=ii*energy*0.1/200.0 + (energy*0.95);
        yy2.push_back(hy->GetSum());
        cout << ii << " " << yy1[ii] << " " << yy2[ii]<< endl;    
    }
    TH1F *peak = new TH1F("@40mI", "Intensity vs. Energy.", 201);
    return peak;
}

void overview_plot(string file, double emin=9000, Double_t emax=11000){
    TCanvas *c1;

    c1=new TCanvas();
    c1->Divide(2,2);
    c1->cd(1);


}
