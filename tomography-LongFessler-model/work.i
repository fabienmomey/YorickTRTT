func lf2d_test2(dir)
{
  if (is_void(dir)) dir = "./";
  db = yhd_restore(dir + "data/DataSet1.yhd");
  p = lf2d_parameters(pixelSize=0.1,
                      pixelStep=0.1,
                      voxelSize=0.1,
                      sourceDistance = 54.1,
                      detectorDistance = 40.8,
                      objectOffset1 = 0.0,
                      objectOffset2 = 0.0,
                      objectDimension1 = 256,
                      objectDimension2 = 256,
                      detectorDimension = 512);
  return lf2d_new_cost_function(param=p, data=db.data, beta=db.theta - pi/2);
}

func lf2d_fix_db(db)
{
  h_set, db,
    pixelSize=0.1,
    pixelStep=0.1,
    voxelSize=0.1,
    sourceDistance = 54.1,
    detectorDistance = 40.8,
    objectOffset = [0.0, 0.0],
    beta = db.theta - pi/2;

}

