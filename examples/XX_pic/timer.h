#ifndef pic_timer_h
#define pic_timer_h

#include <map>
#include <chrono>

// NOTE: This doesn't support pausing
class elapsed_timer_t {
    using time_t = std::chrono::time_point<std::chrono::system_clock>;

    private:
        std::map<std::string,time_t> times;
        std::map<std::string, std::chrono::duration<double> > total_times;

    public:
        void start(std::string label)
        {
            time_t start = std::chrono::system_clock::now();
            times[label] = start;
        }

        double stop(std::string label)
        {
            time_t start = times[label];
            time_t end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end-start;

            std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

            total_times[label] += elapsed_seconds;

            return elapsed_seconds.count();
        }

        double get_total(std::string label)
        {
            return total_times[label].count();
        }
};

class particle_timer_t {

    private:
        const std::string label = "PARTICLES";

        elapsed_timer_t* timer;
    public:

        particle_timer_t() {
            timer = new elapsed_timer_t();
        }

        ~particle_timer_t()
        {
            delete timer;
        }

        void start()
        {
            timer->start(label);
        }

        double stop()
        {
            return timer->stop(label);
        }

        double total()
        {
            return timer->get_total(label);
        }

        double per_particle(size_t count)
        {
            return total() / count;
        }
};

#endif // guard
