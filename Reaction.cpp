#include "Reaction.h"

/*
 * particle_number --
 *   在 particle[] 中查找名为 p 的粒子，若找到则返回索引，否则报错并退出
 */
int particle_number(char* p, char** particle, int n_particles)
{
    for (int i = 0; i < n_particles; i++) {
        // strcmp 返回 0 表示字符串相等
        if (std::strcmp(p, particle[i]) == 0) {
            return i;
        }
    }
    // 未找到
    std::fprintf(stderr, "Error: Particle \"%s\" does not exist.\n", p);
    std::exit(EXIT_FAILURE);
}


/*
 * Read_reaction --
 *   从文件 fp 中读取反应信息，每条反应包含:
 *       反应编号 Rnum[i]
 *       左边粒子 NUM_L 个
 *       右边粒子 NUM_R 个
 *       A, B, E, ER 等
 */
int Read_reaction(
    FILE* fp,
    int NUM_L,
    int NUM_R,
    int n_particles,
    int* Rnum,
    int** reactl,
    int** reactr,
    double* A,
    double* B,
    double* E,
    double* ER,
    char** particle
)
{
    // 先读取行数 n: 跳过 '#' 或空白行
    char line[1000];
    int n = 0;

    // 用 fgets 逐行试读
    while (std::fgets(line, 1000, fp) != NULL) {
        // 若该行是注释行或只有空白字符，则跳过
        if (line[0] == '#' || std::isspace(line[0])) {
            // 保持 n 不变
            continue;
        }
        // 否则 n++
        n++;
    }
    // 回到文件开头
    std::rewind(fp);

    // 逐行读取数据
    int i;
    for (i = 0; i < n; i++) {
        if (std::fgets(line, 1000, fp) == NULL) break;

        // 若行首是 '#' 或空白行，需要跳过
        if (line[0] == '#' || std::isspace(line[0])) {
            // 调整 i-1，以便下一轮循环覆盖这行
            i--;
            continue;
        }

        // 回退文件指针到这一行开头，以用 fscanf 来读
        std::fseek(fp, -static_cast<long>(std::strlen(line)), SEEK_CUR);

        // 读 反应编号
        std::fscanf(fp, "%d", &Rnum[i]);
        // 读 左边 NUM_L 个粒子
        for (int j = 0; j < NUM_L; j++){
            char p[100];
            std::fscanf(fp, "%s", p);
            reactl[Rnum[i]][j] = particle_number(p, particle, n_particles);
        }
        // 读 右边 NUM_R 个粒子
        for (int j = 0; j < NUM_R; j++){
            char p[100];
            std::fscanf(fp, "%s", p);
            reactr[Rnum[i]][j] = particle_number(p, particle, n_particles);
        }

        // 读 4 个数: A,R,B,E
        std::fscanf(fp, "%lf %lf %lf %lf\n",
                    &A[Rnum[i]], &B[Rnum[i]],
                    &E[Rnum[i]], &ER[Rnum[i]]);
    }
    return n;
}


/*
 * Initial_condition --
 *   从 fp 中读取粒子名称和初始浓度，存到 particle[], y[]。
 *   particle[0] = "-" (标记?),  y[0] = 1.0
 *   particle[1] = "M",         y[1] = 0.0
 */
int Initial_condition(FILE* fp, char** particle, double* y)
{
    // 先统计行数
    char line[1000];
    int n = 2;  // 初始多算了2行，以便从第三行开始放真正数据?

    while (std::fgets(line, 1000, fp) != NULL) {
        if (line[0] == '#' || std::isspace(line[0])) {
            continue;
        }
        n++;
    }
    std::rewind(fp);

    // 强行给 0 和 1 做初始值
    std::strcpy(particle[0], "-");
    y[0] = 1.0;
    std::strcpy(particle[1], "M");
    y[1] = 0.0;

    // 从第 i=2 开始往后放
    for (int i = 2; i < n; i++) {
        if (std::fgets(line, 1000, fp) == NULL) break;

        if (line[0] == '#' || std::isspace(line[0])) {
            i--;
            continue;
        }

        int dumn;
        // line里有三项: 整数 dumn, 字符串, 浓度
        std::sscanf(line, "%d %s %lf", &dumn, particle[i], &y[i]);
    }

    return n;
}
