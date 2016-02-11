from gromacs.fileformats import XVG

actual = XVG("actual.xvg")
desired = XVG("desired.xvg")

print(actual.to_df(), desired.to_df())
assert actual.to_df().to_dict() == desired.to_df().to_dict(), actual.to_df() == desired.to_df()
print("mdrun-rerun energies are identical")
