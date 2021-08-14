unit SK42WGS;

interface

uses Math;

procedure WGSToSK42(B, L, H:Double; var X, Y, _H :Double); overload;
procedure WGSToSK42(B, L, H:Double; Zone:Integer; var X, Y, _H :Double); overload;
procedure SK42ToWGS(X, Y, H:Double; var B, L, _H:Double);

implementation

/// Перевод из x, y (Гаусса-Крюгера) в Широту-Долготу (единая СК)

procedure GaussKrugerToGeo_Kras(x, y: double; var OutputB, OutputL: double);

  function XToB(X,Y:double):Double ;
  var
    No : Integer ;
    Bi, Bo, Zo, Ba, Bb, Bc, Db : Double ;
  begin
    No := trunc(Y / 1000000);
    Bi := X / 6367558.4968;

    Bo := Bi + sin(Bi * 2) * (0.00252588685 - 0.0000149186 * power(sin(Bi) , 2)
          + 0.00000011904 * power(sin(Bi) , 4));

    Zo := (Y - (10 * No + 5) * 100000.0) / (6378245.0 * cos(Bo));

    Ba := Zo * Zo * (0.01672 - 0.0063 * power(sin(Bo) , 2)
          + 0.01188 * power(sin(Bo) , 4) - 0.00328 * power(sin(Bo) , 6));
    Bb := Zo * Zo * (0.042858 - 0.025318 * power(sin(Bo) , 2)
          + 0.014346 * power(sin(Bo) , 4) - 0.001264 * power(sin(Bo) , 6) - Ba);
    Bc := Zo * Zo * (0.10500614 - 0.04559916 * power(sin(Bo) , 2)
          + 0.00228901 * power(sin(Bo) , 4)
          - 0.00002987 * power(sin(Bo) , 6) - Bb);

    dB := Zo * Zo * sin(Bo * 2) * (0.251684631
          - 0.003369263 * power(sin(Bo) , 2)
          + 0.000011276 * power(sin(Bo) , 4) - Bc);

    XToB := (Bo - dB) * 180 / Pi;
  end ;

  function YToL(X,Y:double):Double ;
   var No : Integer ;
       Bi, Bo, Zo, La, Lb, Lc, Ld, dL : Double ;
  begin
    No := trunc(Y / 1000000);
    Bi := X / 6367558.4968;
    Bo := Bi + sin(Bi * 2) * (0.00252588685 - 0.0000149186 * power(sin(Bi) , 2)
          + 0.00000011904 * power(sin(Bi) , 4));
    Zo := (Y - (10 * No + 5) * 100000.0) / (6378245.0 * cos(Bo));
    La := Zo * Zo * (0.0038 + 0.0524 * power(sin(Bo) , 2)
          + 0.0482 * power(sin(Bo) , 4) + 0.0032 * power(sin(Bo) , 6));
    Lb := Zo * Zo * (0.01225 + 0.09477 * power(sin(Bo) , 2)
          + 0.03282 * power(sin(Bo) , 4) - 0.00034 * power(sin(Bo) , 6) - La);
    Lc := Zo * Zo * (0.0420025 + 0.1487407 * power(sin(Bo) , 2)
          + 0.005942 * power(sin(Bo) , 4)
          - 0.000015 * power(sin(Bo) , 6) - Lb);
    Ld := Zo * Zo * (0.16778975 + 0.16273586 * power(sin(Bo) , 2)
          - 0.0005249 * power(sin(Bo) , 4)
          - 0.00000846 * power(sin(Bo) , 6) - Lc);
    dL := Zo * (1 - 0.0033467108 * power(sin(Bo) , 2)
          - 0.0000056002 * power(sin(Bo) , 4)
          - 0.0000000187 * power(sin(Bo) , 6) - Ld);
    YToL := (6 * (No - 0.5) / 57.29577951 + dL) * 180 / Pi;
  end;

begin
  OutputB := XToB(x,y);
  OutputL := YToL(x,y);
end;

/// Перевод из Широты-Долготы в x, y (Гаусса-Крюгера) (единая СК)
procedure GeoToGaussKruger_Kras(B1, L1 : double; var Outputx, Outputy: double;
                            var Zone: integer; AutoZone: boolean);   // № зоны, автоматическое определение № зоны вкл/выкл

  procedure GetZone (L :Double);
  begin
    if AutoZone then
      Zone := trunc (6 + L) div  6;
  end;

  function BToX(B, L: DOUBLE ) : double ;
  var No : Integer ;
      Lo, Bo, Xa,Xb,Xc,Xd : double ;
  begin
     No := Zone;
     Lo := (L - (3 + 6 * (No - 1))) / 57.29577951;
     Bo := B * Pi / 180;
     Xa := power(Lo , 2) * (109500.0 - 574700.0 * power(sin(Bo) , 2)
           + 863700.0 * power(sin(Bo) , 4) - 398600.0 * power(sin(Bo) , 6));
     Xb := power(Lo , 2) * (278194.0 - 830174.0 * power(sin(Bo) , 2)
           + 572434.0 * power(sin(Bo) , 4) - 16010.0 * power(sin(Bo) , 6) + Xa);
     Xc := power(Lo , 2) * (672483.4 - 811219.9 * power(sin(Bo) , 2)
           + 5420.0 * power(sin(Bo) , 4) - 10.6 * power(sin(Bo) , 6) + Xb);
     Xd := power(Lo , 2) * (1594561.25 + 5336.535 * power(sin(Bo) , 2)
           + 26.79 * power(sin(Bo) , 4) + 0.149 * power(sin(Bo) , 6) + Xc);
     BToX := 6367558.4968 * Bo - sin(Bo * 2) * (16002.89
           + 66.9607 * power(sin(Bo) , 2) + 0.3515 * power(sin(Bo) , 4) - Xd);
  end;

  function LToY( B,  L: double ): double ;
  var No : Integer;
      Lo, Bo, Ya, Yb, Yc : double ;
  begin
   No := Zone;

   Lo := (L - (3 + 6 * (No - 1))) / 57.29577951;
   Bo := B * Pi / 180;
   Ya := power(Lo , 2) * (79690.0 - 866190.0 * power(sin(Bo) , 2)
       + 1730360.0 * power(sin(Bo) , 4) - 945460.0 * power(sin(Bo) , 6));
   Yb := power(Lo , 2) * (270806.0 - 1523417.0 * power(sin(Bo) , 2)
       + 1327645.0 * power(sin(Bo) , 4) - 21701.0 * power(sin(Bo) , 6) + Ya);
   Yc := power(Lo , 2) * (1070204.16 - 2136826.66 * power(sin(Bo) , 2)
       + 17.98 * power(sin(Bo) , 4) - 11.99 * power(sin(Bo) , 6) + Yb);
   LToY := (5 + 10 * No) * 100000.0 + Lo * cos(Bo) *(6378245
           + 21346.1415 * power(sin(Bo) , 2)
           + 107.159 *power(sin(Bo) , 4) + 0.5977 * power(sin(Bo) , 6) + Yc);
  end ;

begin
  GetZone(L1);
  Outputx := BToX(B1, L1);
  Outputy := LToY(B1, L1);
end;

////////////// перевод Широты-Долготы из одной СК в другую WGS и СК-42, ГОСТ
procedure Geo1ToGeo2(B1, L1, H1 : double;  // Входные широта-долгота-высота
			 Reverse: boolean;
			 var B2, L2, H2 : double);

// Вспомогательные значения для преобразования эллипсоидов

  var
  dx, dy, dz: Double;  // Линейные элементы трансформирования, м
  wx, wy, wz: Double;  // Угловые элементы трансформирования, в секундах
  ms : Double;         // Масштабный коэффициент ppm

  a_1, al_1 : Double;   // Полуось и сжатие  СК 1
  a_2, al_2 : Double;   // Полуось и сжатие  СК 2

  e2_1 : Double;  // Квадрат эксцентриситета эллипсоида 1й СК
  e2_2 : Double;  // Квадрат эксцентриситета эллипсоида 2й СК

  /// Среднее/разности
  a : Double;
  e2 : Double;
  da : Double;
  de2 : Double;

  B2_0, L2_0, H2_0: Double; /// первая итеррация

  const
    Pi = 3.14159265358979;     ///  Число Пи
    ro = 206264.8062471;       ///  Число угловых секунд в радиане

  procedure GetAllConsts;
  begin
    a_1 := 6378245;   al_1 := 1/298.3;
    a_2 := 6378137;   al_2 := 1/298.257223563;

    if Reverse then
    Begin
         a_1 := 6378137;    al_1 := 1/298.257223563;
         a_2 := 6378245;    al_2 := 1/298.3;
    End;

    e2_1 := 2 * al_1 - al_1*al_1;
    e2_2 := 2 * al_2 - al_2*al_2;

    a := (a_1 + a_2) / 2 ;
    e2 := (e2_1 + e2_2) / 2;
    da  := a_1 - a_2  ;
    de2 := e2_1 - e2_2;

	  dx := -23.92;
    dy := 141.27;
    dz := 80.9;

    wx := 0;
    wy := 0.35;
    wz := 0.82;

    ms := 0.12*1E-6;

    if Reverse then
    Begin
         dx := - dx;        dy := - dy;        dz := - dz;
         wx := - wx;        wy := - wy;        wz := - wz;
         ms := - ms;
         a_1 := 6378137;    al_1 := 1/298.257223563;
         a_2 := 6378245;    al_2 := 1/298.3;
    End;
  end;

  function dB(Bd, Ld, H: Double) : Double;
  var B, L, M, N : Double;
  begin
    B := Bd * Pi / 180 ;
    L := Ld * Pi / 180 ;
    M := a * (1 - e2) / Power( (1 - e2 * Sin(B) * Sin(B)), 1.5);
    N := a * Power((1 - e2 * Sin(B) * Sin(B)), -0.5) ;

    dB := ro / (M + H) * (N / a * e2 * Sin(B) * Cos(B) * da
          + (N * N/ (a * a) + 1) * N * Sin(B) * Cos(B) * de2 / 2
          - (dx * Cos(L) + dy * Sin(L)) * Sin(B) + dz * Cos(B))
          - wx * Sin(L) * (1 + e2 * Cos(2 * B))
          + wy * Cos(L) * (1 + e2 * Cos(2 * B))
          - ro * ms * e2 * Sin(B) * Cos(B) ;
  end;

  function dL(Bd, Ld, H: Double) : Double;
  var B, L, N : Double;
  begin
    B := Bd * Pi / 180;
    L := Ld * Pi / 180;
    N := a * Power((1 - e2 * Sin(B) * Sin (B)), -0.5);
    dL := ro / ((N + H) * Cos(B)) * (-dx * Sin(L) + dy * Cos(L))
		      + Tan(B) * (1 - e2) * (wx * Cos(L) + wy * Sin(L)) - wz;
  end;

  function dH(Bd, Ld, H :Double) : Double;
  var B, L, N : Double ;
  Begin
    B := Bd * Pi / 180 ;
    L := Ld * Pi / 180 ;
    N := a * Power((1 - e2 * Sin(B) * Sin(B) ), -0.5) ;
    dH := -a / N * da + N * Sin(B) * Sin(B) * de2 / 2
			    + (dx * Cos(L) + dy * Sin(L)) * Cos(B) + dz * Sin(B)
			    - N * e2 * Sin(B) * Cos(B) * (wx / ro * Sin(L) - wy / ro * Cos(L))
			    + (a * a / N + H) * ms;
  End;

begin

  GetAllConsts;

  B2_0 := dB(B1,L1,H1)/3600;
  L2_0 := dL(B1,L1,H1)/3600;
  H2_0 := dH(B1,L1,H1);

  B2 := dB(B1 - B2_0,L1 - L2_0,H1 - H2_0)/3600;
  L2 := dL(B1 - B2_0,L1 - L2_0,H1 - H2_0)/3600;
  H2 := dH(B1 - B2_0,L1 - L2_0,H1 - H2_0);

  B2 := B1 - (B2_0 + B2)/2;
  L2 := L1 - (L2_0 + L2)/2;
  H2 := H1 - (H2_0 + H2)/2;
end;


procedure WGSToSK42(B, L, H:Double; var X, Y, _H :Double); overload;
var B2, L2, H2:Double;
    zone :Integer;
begin
  Zone := 0;
  Geo1ToGeo2(B, L, H, true, B2, L2, H2);      _H := H2;
  GeoToGaussKruger_Kras(B2, L2, x, y, zone, true);
end;

procedure WGSToSK42(B, L, H:Double; Zone:Integer; var X, Y, _H :Double); overload;
var B2, L2, H2:Double;
begin
  Geo1ToGeo2(B, L, H, true, B2, L2, H2);     _H := H2;
  GeoToGaussKruger_Kras(B2, L2, x, y, zone, false);
end;

procedure SK42ToWGS(X, Y, H:Double; var B, L, _H:Double);
var B2, L2, H2:Double;
begin
  H2 := H;
  GaussKrugerToGeo_Kras(x, y, B2, L2, );
  Geo1ToGeo2(B2, L2, H2, false, B, L, _H);
end;

end.
