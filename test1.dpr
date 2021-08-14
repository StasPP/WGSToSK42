program test1;

{$APPTYPE CONSOLE}

uses
  SysUtils,
  SK42WGS; {������� ���������� ��� ���������� ��42 <-> WGS}
  

var B,L,H, x, y, h1, x2, y2, h2: Double;
begin
  /// ������ �������������

  {�� ��������� ��-42, �������� ������-������� ===> � ������-������� WGS}
  SK42ToWGS(6099100.961, 14628018.567, 175.978,  B, L, H);
            {�����}      {������}      {������}  {���������� ��� ������}

  {�� ������-������� WGS ===> � ��-42, �������� ������-�������}
  WGSToSK42(55.0000, 83.0000, 140.00, x, y, h1);
            {������}{�������}{������}  {���������� ��� ������}

  {�� ������-������� WGS ===> � ��-42, �������� ������-�������, �������� ����}
  WGSToSK42(55.0000, 83.0000, 140.00, 15,     x2, y2, h2);
            {������}{�������}{������}{#����} {���������� ��� ������}

  writeln(formatfloat('00.0000000',B));
  writeln(formatfloat('00.0000000',L));
  writeln(formatfloat('00.000',H));
  writeln;
  writeln(formatfloat('00.000',x));
  writeln(formatfloat('00.000',y));
  writeln(formatfloat('00.000',H1));
  writeln;
  writeln(formatfloat('00.000',x2));
  writeln(formatfloat('00.000',y2));
  writeln(formatfloat('00.000',H2));


  readln;
end.
