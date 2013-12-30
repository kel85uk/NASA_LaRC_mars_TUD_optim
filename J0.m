function j0 = J0(year, month, day)
j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
+ fix(275*month/9) + day + 1721013.5;
end %J0
