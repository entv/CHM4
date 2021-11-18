#include <coroutine>
#include <optional>
#include "LinearAlgebra.h"

class Grid
{
public:
	struct promise_type
	{
		Grid get_return_object()
		{
			return Grid{ Handle::from_promise(*this) };
		}

		static std::suspend_always initial_suspend() noexcept
		{
			return {};
		}

		static std::suspend_always final_suspend() noexcept
		{
			return {};
		}

		std::suspend_always yield_value(Vector point) noexcept
		{
			current_point = std::move(point);
			return {};
		}

		void await_transform() = delete;

		void return_void() noexcept {}

		[[noreturn]]
		static void unhandled_exception()
		{
			throw;
		}

		std::optional<Vector> current_point;
	};

	using Handle = std::coroutine_handle<promise_type>;

	explicit Grid(const Handle coroutine) : m_coroutine{ coroutine }
	{}

	Grid() = default;
	~Grid()
	{
		if (m_coroutine)
		{
			m_coroutine.destroy();
		}
	}

	class Iter
	{
	public:
		void operator++()
		{
			m_coroutine.resume();
		}

		const Vector& operator*() const
		{
			return *m_coroutine.promise().current_point;
		}

		bool operator==(std::default_sentinel_t) const
		{
			return !m_coroutine || m_coroutine.done();
		}

		explicit Iter(const Handle coroutine) :
			m_coroutine{ coroutine }
		{}

	private:
		Handle m_coroutine;
	};

	Iter begin()
	{
		if (m_coroutine)
		{
			m_coroutine.resume();
		}
		return Iter{ m_coroutine };
	}

	std::default_sentinel_t end()
	{
		return {};
	}

private:
	Handle m_coroutine;
};