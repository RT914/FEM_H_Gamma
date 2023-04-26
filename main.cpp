#include <GL/freeglut.h>

#define _USE_MATH_DEFINES
#include <math.h>

//#include "fem.h"
#include "draw.h"
// #include "Camera.h"

int WindowPositionX = 100;  //��������E�B���h�E�ʒu��X���W
int WindowPositionY = 100;  //��������E�B���h�E�ʒu��Y���W
int WindowWidth = 768;    //��������E�B���h�E�̕�
int WindowHeight = 768;    //��������E�B���h�E�̍���
char WindowTitle[] = "FEM";  //�E�B���h�E�̃^�C�g��
int SimulationTime = 0;
bool input_key = false;
int mx, my;

// Camera g_Camera;

//----------------------------------------------------
// �֐��v���g�^�C�v�i��ɌĂяo���֐����ƈ����̐錾�j
//----------------------------------------------------
void Initialize(void);
void Display(void);
void Idle();
void Ground(void);
// void projection_and_modelview(const Camera& in_Camera);
// void mouseDrag(int x, int y);
void mouseDown(int x, int y);
void mouse(int button, int state, int x, int y);
void Keyboard(unsigned char key, int x, int y);

//----------------------------------------------------
// ���C���֐�
//----------------------------------------------------
int main(int argc, char* argv[]) {
	// g_Camera.setEyePoint(Eigen::Vector3d{ 0.0, -100.0, 20.0 });
	// g_Camera.lookAt(Eigen::Vector3d{ 0.0, 500.0, 0.0 }, Eigen::Vector3d{ 0.0, 0.0, 1.0 });

	glutInit(&argc, argv);//���̏�����
	glutInitWindowPosition(WindowPositionX, WindowPositionY); //�E�B���h�E�̈ʒu�̎w��
	glutInitWindowSize(WindowWidth, WindowHeight); //�E�B���h�E�T�C�Y�̎w��
	glutInitDisplayMode(GLUT_RGBA); //�f�B�X�v���C���[�h�̎w��
	glutCreateWindow(WindowTitle);  //�E�B���h�E�̍쐬
	glutIdleFunc(Idle); //�v���O�����A�C�h����Ԏ�(�ɂȎ�)�ɌĂяo�����֐�
	glutKeyboardFunc(Keyboard);//�L�[�{�[�h���͎��ɌĂяo�����֐����w�肷��
	glutDisplayFunc(Display); //�`�掞�ɌĂяo�����֐����w�肷��
	glutMouseFunc(mouse);
	// glutMotionFunc(mouseDrag);
	Initialize(); //�����ݒ�̊֐����Ăяo��
	glutMainLoop();
	return 0;
}
//----------------------------------------------------
// �����ݒ�̊֐�
//----------------------------------------------------
void Initialize(void) {
	glClearColor(1.0, 1.0, 1.0, 1.0); //�w�i�F
	glEnable(GL_DEPTH_TEST);//�f�v�X�o�b�t�@���g�p�FglutInitDisplayMode() �� GLUT_DEPTH ���w�肷��

	gluPerspective(30.0, (double)WindowWidth / (double)WindowHeight, 0.1, 1000.0); //�������e�@�̎��̐�gluPerspactive(th, w/h, near, far);

	//���_�ݒ�
	//gluLookAt(
	//	0.0, -100.0, 20.0, // ���_�̈ʒux,y,z;
	//	0.0, 500.0, 0.0,   // ���E�̒��S�ʒu�̎Q�Ɠ_���Wx,y,z
	//	0.0, 0.0, 1.0);  //���E�̏�����̃x�N�g��x,y,z
}

//----------------------------------------------------
// �`��̊֐�
//----------------------------------------------------
void Display(void) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //�o�b�t�@�̏���
	// projection_and_modelview(g_Camera);

	glColor4f(0.5f, 0.0f, 0.0f, 1.0f);
	//printf("Simulation Time : %d���\n", SimulationTime);
	SimulationTime++;

	//fem(SimulationTime);
	//fem_for_key(SimulationTime, input_key);
	input_key = false;

	glFlush();
};

void Idle() {
	glutPostRedisplay(); //glutDisplayFunc()���P����s����
}

/*
void projection_and_modelview(const Camera& in_Camera)
{
	const double fovy_deg = (2.0 * 180.0 / M_PI) * atan(0.024 * 0.5 / in_Camera.getFocalLength());

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy_deg, double(WindowWidth) / double(WindowHeight), 0.01 * in_Camera.getFocalLength(), 1000.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	const Eigen::Vector3d lookAtPoint = in_Camera.getLookAtPoint();
	gluLookAt(in_Camera.getEyePoint().x(), in_Camera.getEyePoint().y(), in_Camera.getEyePoint().z(), lookAtPoint.x(), lookAtPoint.y(), lookAtPoint.z(), in_Camera.getYVector().x(), in_Camera.getYVector().y(), in_Camera.getYVector().z());
}
*/

/*
void mouseDrag(int x, int y)
{
	const int _dx = x - mx;
	mx = x; my = y;

	const double dx = double(_dx) / double(WindowWidth);
	const double scale = 2.0;

	g_Camera.rotateCameraInLocalFrameFixLookAt(dx * scale);
	glutPostRedisplay();
}
*/

void mouseDown(int x, int y)
{
	mx = x; my = y;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
		mouseDown(x, y);
}

void Keyboard(unsigned char key, int x, int y) {
	switch (key)
	{
	case 'w':
		input_key = true;
	default:
		break;
	}
}