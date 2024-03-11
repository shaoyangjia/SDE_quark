/*
Definition of the class GH51 for the roots and weights of Gauss-Hermite
quadrature.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark.

SDE_quark is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SDE_quark is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SDE_quark. If not, see <https://www.gnu.org/licenses/>.
*/
class GH51
{
public:
static unsigned n_grid;// number of grid points
static double r[];// roots
static double w[];// weights
};

unsigned GH51::n_grid = 51;
double GH51::r[] =
{-9.284352965094795351319589826744050e+00,
-8.627089729363683900942305626813322e+00,
-8.082012673504012312264421780128032e+00,
-7.594737944583593858283165900502354e+00,
-7.144520571479617387922189664095640e+00,
-6.720677095868753170293530274648219e+00,
-6.316797657903404861201579478802159e+00,
-5.928653885873137241446784173604101e+00,
-5.553268005264917483998488023644313e+00,
-5.188439342938137244232166267465800e+00,
-4.832479614578074844644106633495539e+00,
-4.484054079057350428172412648564205e+00,
-4.142080801962062075460835330886766e+00,
-3.805663856751933415978328412165865e+00,
-3.474047361655858257023510304861702e+00,
-3.146582843771695614520922390511259e+00,
-2.822705423305429395242072132532485e+00,
-2.501916004057736842014492140151560e+00,
-2.183767652536192649392887688009068e+00,
-1.867854955699560282056381765869446e+00,
-1.553805529455421829965189317590557e+00,
-1.241273096623611404965004112455063e+00,
-9.299317156001758455374783807201311e-01,
-6.194708497471526076338932398357429e-01,
-3.095910409372403804262319226836553e-01,
0.000000000000000000000000000000000e+00,
3.095910409372403804262319226836553e-01,
6.194708497471526076338932398357429e-01,
9.299317156001758455374783807201311e-01,
1.241273096623611404965004112455063e+00,
1.553805529455421829965189317590557e+00,
1.867854955699560282056381765869446e+00,
2.183767652536192649392887688009068e+00,
2.501916004057736842014492140151560e+00,
2.822705423305429395242072132532485e+00,
3.146582843771695614520922390511259e+00,
3.474047361655858257023510304861702e+00,
3.805663856751933415978328412165865e+00,
4.142080801962062075460835330886766e+00,
4.484054079057350428172412648564205e+00,
4.832479614578074844644106633495539e+00,
5.188439342938137244232166267465800e+00,
5.553268005264917483998488023644313e+00,
5.928653885873137241446784173604101e+00,
6.316797657903404861201579478802159e+00,
6.720677095868753170293530274648219e+00,
7.144520571479617387922189664095640e+00,
7.594737944583593858283165900502354e+00,
8.082012673504012312264421780128032e+00,
8.627089729363683900942305626813322e+00,
9.284352965094795351319589826744050e+00};
double GH51::w[] =
{7.586846833851778049151448612974491e-01,
5.864287425038513523745109523588326e-01,
5.113148070267715672443387120438274e-01,
4.664299739336288941338182212348329e-01,
4.357071372202365178161187486693962e-01,
4.130175315947753333922776164399693e-01,
3.954318931792801250502122911711922e-01,
3.813445580868277118291587157727918e-01,
3.697897868419675604378937805449823e-01,
3.601463743888477786470048158662394e-01,
3.519939135068397839667397875018651e-01,
3.450363958905212635741577287262771e-01,
3.390587302852089979232630412298022e-01,
3.339006067526961718527900302433409e-01,
3.294400717697356650859319415758364e-01,
3.255828179885981166208352988178376e-01,
3.222549814252920197255036782735260e-01,
3.193981705344928201384391286410391e-01,
3.169659611758387507407519478874747e-01,
3.149213821221034659281201584235532e-01,
3.132350877650721399092503816063982e-01,
3.118840198763703774886835162760690e-01,
3.108504266084930467428648626082577e-01,
3.101211499903285817580922412162181e-01,
3.096871220156214854490883681137348e-01,
3.095430294466409715248289558076067e-01,
3.096871220156214854490883681137348e-01,
3.101211499903285817580922412162181e-01,
3.108504266084930467428648626082577e-01,
3.118840198763703774886835162760690e-01,
3.132350877650721399092503816063982e-01,
3.149213821221034659281201584235532e-01,
3.169659611758387507407519478874747e-01,
3.193981705344928201384391286410391e-01,
3.222549814252920197255036782735260e-01,
3.255828179885981166208352988178376e-01,
3.294400717697356650859319415758364e-01,
3.339006067526961718527900302433409e-01,
3.390587302852089979232630412298022e-01,
3.450363958905212635741577287262771e-01,
3.519939135068397839667397875018651e-01,
3.601463743888477786470048158662394e-01,
3.697897868419675604378937805449823e-01,
3.813445580868277118291587157727918e-01,
3.954318931792801250502122911711922e-01,
4.130175315947753333922776164399693e-01,
4.357071372202365178161187486693962e-01,
4.664299739336288941338182212348329e-01,
5.113148070267715672443387120438274e-01,
5.864287425038513523745109523588326e-01,
7.586846833851778049151448612974491e-01};
