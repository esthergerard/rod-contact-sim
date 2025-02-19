#include "world.h"

world::world()
{
    ;
}

world::world(setInput &m_inputData)
{
    render = m_inputData.GetBoolOpt("render");     // boolean
    saveData = m_inputData.GetBoolOpt("saveData"); // boolean

    // Physical parameters
    helixPitch = m_inputData.GetScalarOpt("helixPitch");   // meter
    helixRadius = m_inputData.GetScalarOpt("helixRadius"); // meter
    gVector = m_inputData.GetVecOpt("gVector");            // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");            // maximum number of iterations
    rodRadius = m_inputData.GetScalarOpt("rodRadius");     // meter
    youngM = m_inputData.GetScalarOpt("youngM");           // Pa
    Poisson = m_inputData.GetScalarOpt("Poisson");         // dimensionless
    deltaTime = m_inputData.GetScalarOpt("deltaTime");     // seconds
    totalTime = m_inputData.GetScalarOpt("totalTime");     // seconds
    tol = m_inputData.GetScalarOpt("tol");                 // small number like 10e-7
    stol = m_inputData.GetScalarOpt("stol");               // small number, e.g. 0.1%
    density = m_inputData.GetScalarOpt("density");         // kg/m^3
    viscosity = m_inputData.GetScalarOpt("viscosity");     // viscosity in Pa-s
    epsilon = m_inputData.GetScalarOpt("epsilon");
    numRod = m_inputData.GetIntOpt("numFlagella");
    distance = m_inputData.GetScalarOpt("distance"); // meter
    omega = m_inputData.GetScalarOpt("omega");       // rad/s

    axisLengthInput = m_inputData.GetScalarOpt("axisLengthInput");
    deltaLengthInput = m_inputData.GetScalarOpt("deltaLengthInput");

    delta = m_inputData.GetScalarOpt("delta");        // meter
    col_limit = m_inputData.GetScalarOpt("colLimit"); // meter
    mu = m_inputData.GetScalarOpt("mu");
    nu = m_inputData.GetScalarOpt("nu"); // m/s
    line_search = m_inputData.GetIntOpt("lineSearch");
    ipc = m_inputData.GetIntOpt("ipc");

    pull_time = 2;    // get time of pulling
    release_time = 0; // get time of loosening
    wait_time = 0;    // get time of waiting
    pull_speed = 0.06;

    ////////////////
    // rod length and numVertices defined in function world::rodGeometry()
    //  double nTurn = axisLengthInput / helixPitch;
    //  RodLength = nTurn * sqrt((2 * M_PI * helixRadius) * (2 * M_PI * helixRadius) + helixPitch * helixPitch);
    //  int newNe = RodLength / deltaLengthInput;
    //  numVertices = newNe + 1;

    ///////////////////

    shearM = youngM / (2.0 * (1.0 + Poisson)); // shear modulus
}

world::~world()
{
    ;
}

bool world::isRender()
{
    return render;
}

void world::OpenFile(ofstream &simfile, ofstream &configfile, string file_type)
{

    int systemRet = system("mkdir -p datafiles"); // make the directory
    if (systemRet == -1)
    {
        cout << "Error in creating directory\n";
    }

    ostringstream simfile_name;
    simfile_name.precision(6);
    simfile_name << "datafiles/simData";
    simfile_name << "_numFlagella_" << numRod;
    simfile_name << "_EA_" << youngM;
    simfile_name << "_dt_" << deltaTime;
    simfile_name << "_totalTime_" << totalTime;
    simfile_name << "_imc_" << (!ipc);
    simfile_name << "_friction_" << mu;
    simfile_name << ".txt";

    simfile.open(simfile_name.str().c_str());
    simfile.precision(10);

    ostringstream configfile_name;
    configfile_name.precision(6);
    configfile_name << "datafiles/configData";
    configfile_name << "_numFlagella_" << numRod;
    configfile_name << "_EA_" << youngM;
    configfile_name << "_dt_" << deltaTime;
    configfile_name << "_totalTime_" << totalTime;
    configfile_name << "_imc_" << (!ipc);
    configfile_name << "_friction_" << mu;
    configfile_name << ".txt";

    configfile.open(configfile_name.str().c_str());
    configfile.precision(10);
}

void world::CloseFile(ofstream &simfile, ofstream &configfile)
{
    if (saveData == false)
        return;
    simfile.close();
    configfile.close();
}

void world::CoutData(ofstream &simfile, ofstream &configfile, double &time_taken)
{
    if (saveData == false)
    {
        return;
    }
    double contact_stiffness;
    if (!ipc)
    {
        contact_stiffness = m_contactPotentialIMC->contact_stiffness;
    }
    else
    {
        contact_stiffness = m_contactPotentialIPC->contact_stiffness;
    }
    // calculate the mean value of the end distance
    double dis = 0;
    int count = 0;
    vector<Vector3d> endNode;
    endNode.clear();
    for (int n = 0; n < rodsVector.size(); n++)
    {
        endNode.push_back(rodsVector[n]->getVertex(rodsVector[n]->nv - 1));
    }

    for (int i = 0; i < endNode.size() - 1; i++)
    {
        Vector3d temp = endNode[i];
        for (int j = i + 1; j < endNode.size(); j++)
        {
            Vector3d temp1 = endNode[j];
            dis = dis + (temp - temp1).norm();
            count++;
        }
    }
    dis = dis / count;

    // external forces (should be equal to the propulisve force)
    Vector3d Fp(0, 0, 0);
    for (int n = 0; n < rodsVector.size(); n++)
    {
        Fp = Fp + stepper->force.segment(rodsVector[n]->ndof * n, 3) + stepper->force.segment(rodsVector[n]->ndof * n + 4, 3);
    }

    // // hydrodynamic force
    // Vector3d fp(0, 0, 0);
    // for (int i = 0; i < numRod * numVertices; i++) {
    //     fp(0) = fp(0) + m_RegularizedStokeslet->ViscousForce(3 * i);
    //     fp(1) = fp(1) + m_RegularizedStokeslet->ViscousForce(3 * i + 1);
    //     fp(2) = fp(2) + m_RegularizedStokeslet->ViscousForce(3 * i + 2);
    // }
    // inertial force
    Vector3d inertiaF(0, 0, 0);
    for (int n = 0; n < numRod; n++)
    {
        for (int i = 0; i < rodsVector[n]->nv; i++)
        {
            inertiaF(0) = inertiaF(0) + v_inertialForce[n]->ForceVec[4 * i];
            inertiaF(1) = inertiaF(1) + v_inertialForce[n]->ForceVec[4 * i + 1];
            inertiaF(2) = inertiaF(2) + v_inertialForce[n]->ForceVec[4 * i + 2];
        }
    }

    simfile << currentTime << " " << iter << " " << m_collisionDetector->num_collisions << " " << time_taken << " " << m_collisionDetector->min_dist << " " << contact_stiffness << " " << dis << " " << Fp(0) << " " << Fp(1) << " " << Fp(2) << " " << inertiaF(0) << " " << inertiaF(1) << " " << inertiaF(2) << " "
            << normf << " "
            << normf0 << endl;

    // Record only every second 0.1 seconds so files don't get humongous
    int rate = int(0.1 / deltaTime);
    if (timeStep % rate == 0)
    {
        for (int n = 0; n < numRod; n++)
        {
            for (int i = 0; i < rodsVector[n]->nv; i++)
            {
                Vector3d xCurrent = rodsVector[n]->getVertex(i);
                configfile << n << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
            }
        }
    }
}

void world::setRodStepper()
{
    // set up geometry
    for (int i = 0; i < numRod; i++)
    {

        // set up each rod
        rodGeometry(i);

        // create the rods
        rodsVector.push_back(new elasticRod(vertices, vertices, density, rodRadius, deltaTime, youngM, shearM, RodLength, friction));

        // set up boundary conditions
        rodBoundaryCondition(i);

        // set up the rod so that all the relevant varibles are populated
        rodsVector[i]->setup();

        // allocate for first iteration
        rodsVector[i]->updateTimeStep();
    }
    // // define the pattern
    // double theta = (numRod - 2) * 1.0 * M_PI / (2 * numRod);
    // double R = distance / 2.0 / cos(theta);
    // theta = 2 * M_PI / numRod;
    // for (int i = 0; i < numRod; i++) {

    //     double t = theta * i;
    //     double offset_x = R - R * cos(t);
    //     double offset_y = R * sin(t);
    //     // start from (0, 0)
    //     rodGeometry(offset_x, offset_y, 0);

    //     rodsVector.push_back(new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
    //                                         youngM, shearM, RodLength, friction));

    //     rodBoundaryCondition(i);
    //     rodsVector[i]->setup();
    //     rodsVector[i]->updateTimeStep();
    // }

    // define the timestepper solver
    stepper = new timeStepper(rodsVector);
    m_collisionDetector = new collisionDetector(rodsVector, *stepper, delta, col_limit);
    if (!ipc)
    {
        m_contactPotentialIMC = new contactPotentialIMC(rodsVector, *stepper, *m_collisionDetector, delta, mu, nu);
    }
    else
    {
        m_contactPotentialIPC = new contactPotentialIPC(rodsVector, *stepper, *m_collisionDetector, delta);
    }

    dx = stepper->dx;

    for (int i = 0; i < numRod; i++)
    {
        v_stretchForce.push_back(new elasticStretchingForce(*rodsVector[i], *stepper, i));
        v_bendingForce.push_back(new elasticBendingForce(*rodsVector[i], *stepper, i));
        v_twistingForce.push_back(new elasticTwistingForce(*rodsVector[i], *stepper, i));
        v_inertialForce.push_back(new inertialForce(*rodsVector[i], *stepper, i));
        v_gravityForce.push_back(new externalGravityForce(*rodsVector[i], *stepper, gVector, i));
    }
    // // define RSS model
    // m_RegularizedStokeslet = new RegularizedStokeslet(rodsVector, *stepper, viscosity, epsilon);

    timeStep = 0;
    currentTime = 0.0;

    Nstep = totalTime / deltaTime;

    // Find out the tolerance, e.g. how small is enough?
    characteristicForce = M_PI * pow(rodRadius, 4) / 4.0 * youngM / pow(RodLength, 2) * numRod;
    forceTol = tol * characteristicForce;

    ne = numVertices - 1;
    nv = numVertices;
}

void world::rodGeometry(int &index_Rod)
{
    // Define the input file containing the knot config and open it
    // string knot_config = "reef" + std::to_string(index_Rod + 1);
    string knot_config = "n1";
    ifstream myfile(("knot_configurations/" + knot_config).c_str());
    if (!myfile.is_open())
    {
        std::cerr << "Failed to open the file." << std::endl;
    }

    // Read the knot config file to extract the vertices values
    int lineCount = 0;
    std::string line;
    std::vector<std::vector<double>> matrix_a;

    while (std::getline(myfile, line))
    {
        double a;
        std::vector<double> vector_a;
        istringstream iss(line);
        while (iss >> a)
        {
            vector_a.push_back(a);
        }
        matrix_a.push_back(vector_a);
        // Increment the line count for each line read
        lineCount++;
    }
    numVertices = lineCount;

    // Create the Eigen matrix structure containing the values of the input file
    int ColNb = matrix_a[0].size();
    MatrixXd data = MatrixXd(numVertices, ColNb);
    for (int RowIdx = 0; RowIdx < numVertices; ++RowIdx)
    {
        for (int ColIdx = 0; ColIdx < ColNb; ++ColIdx)
        {
            data(RowIdx, ColIdx) = matrix_a[RowIdx][ColIdx];
        }
    }

    // INPUT DATA POSITIONS ARE SUPPOSED TO BE IN MM
    vertices = MatrixXd(numVertices, 3);
    if (ColNb == 3)
    {
        // vertices = data * 1e-3;
        vertices = data;
    }
    else
    {
        for (int RowIdx = 0; RowIdx < numVertices; ++RowIdx)
        {
            for (int ColIdx = 0; ColIdx < 3; ++ColIdx)
            {
            //     vertices(RowIdx, ColIdx) = data(RowIdx, ColIdx) * 1e-3;
                vertices(RowIdx, ColIdx) = data(RowIdx, ColIdx);
            }
        }
    }

    // RodLength definition : euclidian distance from the first to the last node coordinates
    RodLength = 0;
    for (int RowIdx = 0; RowIdx < numVertices - 1; ++RowIdx)
    {
        RodLength += std::sqrt((vertices(RowIdx + 1, 0) - vertices(RowIdx, 0)) * (vertices(RowIdx + 1, 0) - vertices(RowIdx, 0)) +
                               (vertices(RowIdx + 1, 1) - vertices(RowIdx, 1)) * (vertices(RowIdx + 1, 1) - vertices(RowIdx, 1)) +
                               (vertices(RowIdx + 1, 2) - vertices(RowIdx, 2)) * (vertices(RowIdx + 1, 2) - vertices(RowIdx, 2)));
    }

    // // arrange the rotation and position of each flagella initially
    // vertices = MatrixXd(numVertices, 3);

    // double helixA = helixRadius;
    // double helixB = helixPitch / (2.0 * M_PI);

    // double T = axisLengthInput / helixB;

    // double newRodLength = T * sqrt(helixA * helixA + helixB * helixB);
    // // define rod length and numVertices
    // RodLength = newRodLength;

    // int newNe = RodLength / deltaLengthInput;
    // double dl = RodLength / newNe;
    // RodLength = RodLength + dl + helixA * sqrt(2);

    // int addVertices = helixA * sqrt(2) / dl;

    // numVertices = newNe + 1 + addVertices + 1;

    // vertices = MatrixXd::Zero(numVertices, 3);
    // vertices(0, 0) = -dl - helixA;
    // vertices(0, 1) = 0 + offset_x;
    // vertices(0, 2) = 0 + offset_y;

    // vertices(1, 0) = -helixA;
    // vertices(1, 1) = 0 + offset_x;
    // vertices(1, 2) = 0 + offset_y;

    // for (int i = 1; i < 2 + addVertices; i++) {
    //     vertices(i, 0) = -helixA + (i - 1) * helixA / addVertices;
    //     vertices(i, 1) = (i - 1) * helixA / addVertices + offset_x;
    //     vertices(i, 2) = 0 + offset_y;
    // }

    // double t = T / newNe;

    // for (int i = 2 + addVertices - 1; i < numVertices; i++) {
    //     T = (i - (2 + addVertices - 1)) * t;
    //     vertices(i, 0) = helixB * T;
    //     vertices(i, 1) = helixA * cos(T) + offset_x;
    //     vertices(i, 2) = helixA * sin(T) + offset_y;
    // }

    // Matrix3d R;
    // R << 1, 0, 0,
    //         0, cos(offset_theta), -sin(offset_theta),
    //         0, sin(offset_theta), cos(offset_theta);

    // for (int i = 0; i < numVertices; i++) {
    //     vertices.row(i) = vertices.row(i) * R.transpose();
    // }
}

void world::rodBoundaryCondition(int n)
{
    // rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(0), 0);
    // rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(1), 1);
    // rodsVector[n]->setThetaBoundaryCondition(rodsVector[n]->getTheta(0), 0);

    rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(0), 0);
    rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(1), 1);
    rodsVector[n]->setThetaBoundaryCondition(0.0, 0);
    rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(numVertices - 1), numVertices - 1);
    rodsVector[n]->setVertexBoundaryCondition(rodsVector[n]->getVertex(numVertices - 2), numVertices - 2);
    rodsVector[n]->setThetaBoundaryCondition(0.0, numVertices - 2);
}

void world::updateBoundary()
{
    // apply omega for b.c.
    for (int i = 0; i < numRod; i++)
    {

        Vector3d u;
        u(0) = 0;
        u(1) = pull_speed;
        u(2) = 0;

        if (currentTime > wait_time && currentTime <= wait_time + pull_time)
        { // Pulling
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(0) - u * deltaTime, 0);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(1) - u * deltaTime, 1);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(numVertices - 1) + u * deltaTime, numVertices - 1);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(numVertices - 2) + u * deltaTime, numVertices - 2);
        }
        else if (currentTime > wait_time + pull_time &&
                 currentTime <= wait_time + pull_time + release_time)
        { // Loosening
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(0) + u * deltaTime, 0);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(1) + u * deltaTime, 1);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(numVertices - 1) - u * deltaTime, numVertices - 1);
            rodsVector[i]->setVertexBoundaryCondition(rodsVector[i]->getVertex(numVertices - 2) - u * deltaTime, numVertices - 2);
        }

        // deltaTheta = omega * (2.0 * M_PI / 60.0) * deltaTime * (-1.0); //CCW direction

        // rodsVector[i]->x = rodsVector[i]->x0;
        // rodsVector[i]->setThetaBoundaryCondition(rodsVector[i]->getTheta(0) + deltaTheta, 0);
        // rodsVector[i]->updateGuess();
    }
}

void world::updateTimeStep()
{
    alpha = 1;
    bool goodSolved = false;
    updateBoundary();
    newtonMethod(goodSolved);

    // calculate pull forces;
    calculateForce();

    // update time step
    for (int i = 0; i < numRod; i++)
    {
        rodsVector[i]->updateTimeStep();
    }
    printSimData();

    currentTime += deltaTime;

    timeStep++;
}

void world::calculateForce() {
    stepper->setZero();
    
    for (int n = 0; n < numRod; n++){
        rodsVector[n]->prepareForIteration();
        v_inertialForce[n]->computeFi();
        v_stretchForce[n]->computeFs();
        v_bendingForce[n]->computeFb();
        v_twistingForce[n]->computeFt();
        v_gravityForce[n]->computeFg();

        temp[0] = stepper->force[0] + stepper->force[4];
        temp[1] = stepper->force[1] + stepper->force[5];
        temp[2] = stepper->force[2] + stepper->force[6];

        temp1[0] = stepper->force[rodsVector[n]->ndof - 3] + stepper->force[rodsVector[n]->ndof - 7];
        temp1[1] = stepper->force[rodsVector[n]->ndof - 2] + stepper->force[rodsVector[n]->ndof - 6];
        temp1[2] = stepper->force[rodsVector[n]->ndof - 1] + stepper->force[rodsVector[n]->ndof - 5];
    }
}

void world::newtonDamper()
{
    if (iter < 10)
        alpha = 1.0;
    else
        alpha *= 0.90;

    if (alpha < 0.1)
        alpha = 0.1;
}

void world::newtonMethod(bool &solved)
{
    normf = forceTol * 10.0;
    normf0 = 0;

    iter = 0;

    // m_RegularizedStokeslet->prepareForViscousForce();

    // Use this only at the start of sim to initialize
    m_collisionDetector->constructCandidateSet();

    while (solved == false)
    {
        stepper->setZero();
        for (int n = 0; n < numRod; n++)
        {
            rodsVector[n]->prepareForIteration();
            v_inertialForce[n]->computeFi();
            v_inertialForce[n]->computeJi();

            v_stretchForce[n]->computeFs();
            v_stretchForce[n]->computeJs();

            v_bendingForce[n]->computeFb();
            v_bendingForce[n]->computeJb();

            v_twistingForce[n]->computeFt();
            v_twistingForce[n]->computeJt();

            v_gravityForce[n]->computeFg();
            v_gravityForce[n]->computeJg();
        }

        // m_RegularizedStokeslet->computeFrs(); // RSS model

        m_collisionDetector->detectCollisions();
        if (!ipc)
        {
            if (iter == 0)
            {
                m_contactPotentialIMC->updateContactStiffness();
            }
            m_contactPotentialIMC->computeFcJc();
        }
        else
        {
            if (iter == 0)
            {
                m_contactPotentialIPC->updateContactStiffness();
            }
            m_contactPotentialIPC->computeFcJc();
        }

        // Compute norm of the force equations.
        normf = stepper->Force.norm();

        if (iter == 0)
            normf0 = normf;

        if (std::isnan(normf) || std::isnan(normf0))
        {
            cout << "Error. NaN value detected. Exiting.\n";
            break;
        }

        if (normf <= forceTol || iter > 0 && normf <= normf0 * stol)
        {
            // if (normf < 1e-2 * m_RegularizedStokeslet->getViscousForceNorm() || normf < sqrt(numRod) * 1e-5) {
            solved = true;
            iter++;
            // }
        }

        if (solved == false)
        {
            stepper->integrator(); // Solve equations of motion

            if (line_search)
            {
                if (!ipc)
                {
                    alpha = lineSearch();
                }
                else
                {
                    alpha = IPCLineSearch();
                }
            }
            else
            {
                newtonDamper();
            }

            for (int n = 0; n < numRod; n++)
            {
                rodsVector[n]->updateNewtonX(dx, n, alpha);
            }
            iter++;
        }

        if (iter > maxIter)
        {
            cout << "Error. Could not converge. Exiting.\n";
            break;
        }
    }

    if (solved == false)
    {
        timeStep = Nstep; // we are exiting
    }
}

void world::printSimData()
{
    double contact_stiffness;
    if (!ipc)
    {
        contact_stiffness = m_contactPotentialIMC->contact_stiffness;
    }
    else
    {
        contact_stiffness = m_contactPotentialIPC->contact_stiffness;
        mu = 0.0;
    }

    printf("time: %.4f | iters: %i | con: %i | min_dist: %.6f | k: %.3e | fric: %.2f\n",
           currentTime, iter,
           m_collisionDetector->num_collisions,
           m_collisionDetector->min_dist,
           contact_stiffness,
           mu);
}

int world::simulationRunning()
{
    if (timeStep < Nstep)
        return 1;
    else
    {
        return -1;
    }
}

int world::numPoints()
{
    return rodsVector[0]->nv;
}

double world::getScaledCoordinate(int i, int j)
{
    return rodsVector[i]->x[j] / RodLength / 2;
}

double world::getCurrentTime()
{
    return currentTime;
}

double world::getTotalTime()
{
    return totalTime;
}

double world::IPCLineSearch()
{
    if (m_collisionDetector->num_collisions == 0)
        return 1;
    double alpha = m_contactPotentialIPC->computeUpperBound(-stepper->DX);
    double alpha_max = 0.95 * alpha;
    // store the current poses
    for (int n = 0; n < numRod; n++)
    {
        rodsVector[n]->xold = rodsVector[n]->x;
    }
    double q0 = stepper->Force.norm();
    alpha = 2 * alpha_max;

    do
    {
        alpha = alpha / 2.0;
        stepper->setZero();
        for (int n = 0; n < numRod; n++)
        {
            rodsVector[n]->x = rodsVector[n]->xold;
            rodsVector[n]->updateNewtonX(dx, n, alpha);

            rodsVector[n]->prepareForIteration();
            v_inertialForce[n]->computeFi();
            v_stretchForce[n]->computeFs();
            v_bendingForce[n]->computeFb();
            v_twistingForce[n]->computeFt();
            v_gravityForce[n]->computeFg();
        }

        // m_RegularizedStokeslet->computeFrs();
        m_collisionDetector->detectCollisions();
        m_contactPotentialIPC->computeFc();

        if (alpha < 1e-5)
            break;
    } while (stepper->Force.norm() >= q0);

    alpha = (alpha < alpha_max) ? alpha : alpha_max;

    for (int n = 0; n < numRod; n++)
    {
        rodsVector[n]->x = rodsVector[n]->xold;
    }

    return alpha;
}

double world::lineSearch()
{
    // store current x
    for (int n = 0; n < numRod; n++)
    {
        rodsVector[n]->xold = rodsVector[n]->x;
    }
    // Initialize an interval for optimal learning rate alpha
    double amax = 2;
    double amin = 1e-3;
    double al = 0;
    double au = 1;

    double a = 1;

    // compute the slope initially
    double q0 = 0.5 * pow(stepper->Force.norm(), 2); // should decrease with the updated learning rate from lineSearch
    double dq0 = -(stepper->Force.transpose() * stepper->Jacobian * stepper->DX)(0);

    bool success = false;
    double m2 = 0.9;
    double m1 = 0.1;
    int iter_l = 0;
    while (!success)
    {

        stepper->setZero();
        for (int n = 0; n < numRod; n++)
        {
            rodsVector[n]->x = rodsVector[n]->xold;
            rodsVector[n]->updateNewtonX(dx, n, a);

            rodsVector[n]->prepareForIteration();
            v_inertialForce[n]->computeFi();
            v_stretchForce[n]->computeFs();
            v_bendingForce[n]->computeFb();
            v_twistingForce[n]->computeFt();
            v_gravityForce[n]->computeFg();
        }

        // m_RegularizedStokeslet->computeFrs();
        m_collisionDetector->detectCollisions();
        m_contactPotentialIMC->computeFc();

        double q = 0.5 * pow(stepper->Force.norm(), 2);
        double slope = (q - q0) / a;

        // if (q != q) {
        //     cout << "cest q" << q << endl;
        // }

        if (slope >= m2 * dq0 && slope <= m1 * dq0)
        {
            success = true;
        }
        else
        {
            if (slope < m2 * dq0)
            {
                al = a;
            }
            else
            {
                au = a;
            }

            if (au < amax)
            {
                a = 0.5 * (al + au);
            }
        }
        if (a > amax || a < amin || iter_l > 100)
        {
            break;
        }
        iter_l++;
    }

    for (int n = 0; n < numRod; n++)
    {
        rodsVector[n]->x = rodsVector[n]->xold;
    }
    return a;
}
