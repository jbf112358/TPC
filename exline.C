void exline()
{
//=========Macro generated from canvas: c4/
//=========  (Sat Dec 28 10:05:59 2019) by ROOT version 6.18/00
   TCanvas *c4 = new TCanvas("c4", "",652,101,600,600);
   gStyle->SetOptFit(1);
   c4->Range(-1.105983,-1.099179,1.105983,1.099179);
   TView *view1 = TView::CreateView(1);
   view1->SetRange(0,0,0,100,100,140);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetTheta(50.34312);
   c4->SetPhi(-80.63664);
   c4->SetLeftMargin(0.2);
   c4->SetRightMargin(0.2);
   c4->SetFrameBorderMode(0);
   
   TGraph2D *graph2d = new TGraph2D(49);
   graph2d->SetName("Graph2D");
   graph2d->SetTitle(";x/mm;y/mm;z/mm");
   graph2d->SetDirectory(0);
   graph2d->GetXaxis()->SetTitleOffset(1.28);
   graph2d->GetXaxis()->SetTitleSize(0.05);
   graph2d->GetXaxis()->SetLabelSize(0.05);
   graph2d->GetXaxis()->CenterTitle(true);
   graph2d->GetXaxis()->SetNdivisions(505);
   graph2d->GetYaxis()->SetTitleOffset(1.28);
   graph2d->GetYaxis()->SetTitleSize(0.05);
   graph2d->GetYaxis()->SetLabelSize(0.05);
   graph2d->GetYaxis()->CenterTitle(true);
   graph2d->GetYaxis()->SetNdivisions(505);
   graph2d->GetZaxis()->SetTitleOffset(1.18);
   graph2d->GetZaxis()->SetTitleSize(0.05);
   graph2d->GetZaxis()->SetLabelSize(0.05);
   graph2d->GetZaxis()->CenterTitle(true);

   graph2d->SetPoint(0,0,0,0);
   graph2d->SetPoint(1,100,100,140);
   graph2d->SetPoint(2,32.49088,100,79.2);
   graph2d->SetPoint(3,31.58367,97.5,78.2);
   graph2d->SetPoint(4,31.40223,95,78);
   graph2d->SetPoint(5,31.03935,92.5,77.6);
   graph2d->SetPoint(6,30.94863,90,77.5);
   graph2d->SetPoint(7,30.31359,87.5,76.8);
   graph2d->SetPoint(8,31.24995,86.67375,77);
   graph2d->SetPoint(9,30.40431,85,76.9);
   graph2d->SetPoint(10,30.13214,82.5,76.6);
   graph2d->SetPoint(11,29.31566,80,75.7);
   graph2d->SetPoint(12,29.16662,79.6068,76);
   graph2d->SetPoint(13,29.58782,77.5,76);
   graph2d->SetPoint(14,29.13422,75,75.5);
   graph2d->SetPoint(15,28.68062,72.5,75);
   graph2d->SetPoint(16,28.5899,70,74.9);
   graph2d->SetPoint(17,28.5899,67.5,74.9);
   graph2d->SetPoint(18,27.08329,65.4729,74);
   graph2d->SetPoint(19,27.77341,65,74);
   graph2d->SetPoint(20,27.59197,62.5,73.8);
   graph2d->SetPoint(21,26.95693,60,73.1);
   graph2d->SetPoint(22,26.68477,57.5,72.8);
   graph2d->SetPoint(23,26.32189,55,72.4);
   graph2d->SetPoint(24,25.59612,52.5,71.6);
   graph2d->SetPoint(25,25.23324,50,71.2);
   graph2d->SetPoint(26,24.99996,48.51222,71.6);
   graph2d->SetPoint(27,25.23324,47.5,71.2);
   graph2d->SetPoint(28,24.96108,45,70.9);
   graph2d->SetPoint(29,24.32604,42.5,70.2);
   graph2d->SetPoint(30,24.1446,40,70);
   graph2d->SetPoint(31,23.50955,37.5,69.3);
   graph2d->SetPoint(32,23.14667,35,68.9);
   graph2d->SetPoint(33,22.87451,32.5,68.6);
   graph2d->SetPoint(34,22.96523,30,68.7);
   graph2d->SetPoint(35,22.91663,28.80866,68.1);
   graph2d->SetPoint(36,22.33019,27.5,68);
   graph2d->SetPoint(37,22.14875,25,67.8);
   graph2d->SetPoint(38,21.33226,22.5,66.9);
   graph2d->SetPoint(39,21.0601,20,66.6);
   graph2d->SetPoint(40,20.69722,17.5,66.2);
   graph2d->SetPoint(41,20.06218,15,65.5);
   graph2d->SetPoint(42,19.97146,12.5,65.4);
   graph2d->SetPoint(43,19.88074,10,65.3);
   graph2d->SetPoint(44,20.8333,9.855596,65.1);
   graph2d->SetPoint(45,18.74997,7.960289,64.8);
   graph2d->SetPoint(46,19.42713,7.5,64.8);
   graph2d->SetPoint(47,18.79209,5,64.1);
   graph2d->SetPoint(48,18.97353,2.5,64.3);
   graph2d->Draw("p0");
   
   TPolyLine3D *pline3D = new TPolyLine3D(2,"");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   pline3D->SetLineColor(ci);
   pline3D->SetLineWidth(3);
   pline3D->SetPoint(0,32.26791,100,78.88782);
   pline3D->SetPoint(1,25.6849,50,71.80824);
   pline3D->Draw();
   
   pline3D = new TPolyLine3D(2,"");

   ci = TColor::GetColor("#0000ff");
   pline3D->SetLineColor(ci);
   pline3D->SetLineWidth(3);
   pline3D->SetPoint(0,25.44855,50,71.51118);
   pline3D->SetPoint(1,18.36107,0,63.51204);
   pline3D->Draw();
   c4->Modified();
   c4->cd();
   c4->SetSelected(c4);
}
